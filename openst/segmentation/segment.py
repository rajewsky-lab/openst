import h5py
import logging
import dask
import dask_image
import dask_image.imread
from dask_image.ndmeasure._utils._label import (
        connected_components_delayed,
        label_adjacency_graph,
        relabel_blocks,
    )
import dask.array as da
from dask.diagnostics import ProgressBar
import os
import pathlib

import numpy as np
from PIL import Image
from skimage import measure
from skimage.io import imsave

from ome_zarr.io import parse_url
from ome_zarr.writer import write_image
import zarr
from openst.utils.file import check_directory_exists, check_file_exists
from openst.utils.pimage import mask_tissue
from openst.utils.pseudoimage import create_unpaired_pseudoimage
from skimage.segmentation import find_boundaries

# Models pretrained by the Rajewsky lab
OPENST_MODEL_NAMES = [
    "HE_cellpose_rajewsky"
]

# adapted from cellpose
MODEL_DIR = pathlib.Path.home().joinpath(".openst", "cellpose", "models")
_MODEL_URL = "http://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/models"

def cache_model_path(basename):
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    url = f"{_MODEL_URL}/{basename}"
    cached_file = os.fspath(MODEL_DIR.joinpath(basename))
    if not os.path.exists(cached_file):
        from cellpose.utils import download_url_to_file
        logging.info('Downloading: "{}" to {}\n'.format(url, cached_file))
        download_url_to_file(url, cached_file, progress=True)
    return cached_file

# TODO: implement gray_dilation
def assemble_from_tiles(tiles, width: int, height: int, channels: int, tile_size: int = 512):
    """
    Assemble an image from a list of tiles.

    Args:
        tiles (list): List of tiles to assemble into an image.
        width (int): Width of the assembled image.
        height (int): Height of the assembled image.
        channels (int): Number of image channels.
        tile_size (int, optional): Size of the tiles. Defaults to 512.

    Returns:
        numpy.ndarray: Assembled image.

    Notes:
        - This function stitches tiles into an image.
    """

    img_assembled = np.zeros((width, height, channels))

    _image_index = 0
    for x in range(0, width - tile_size, tile_size):
        for y in range(0, height - tile_size, tile_size):
            img_assembled[x[0] : (x[0] + tile_size), y[1] : (y[1] + tile_size)] = tiles[_image_index]
            _image_index += 1

    return img_assembled


def da_imsave(fname, arr, compute=True):
    """Write arr to image

    from https://github.com/dask/dask-image/issues/110
    
    Parameters
    ----------
    fnames: string
        A formatting string like 'myfile{:02d}.png'
        Should support arr.ndims-2 indices to be formatted
    arr: dask.array
        Array of at least 2 dimensions to be written to disk as images
    compute: Boolean (optional)
        whether to write to disk immediately or return a dask.array of the to be written indices
    """
    indices = [da.arange(n, chunks=c) for n,c in zip(arr.shape[:-2], arr.chunksize[:-2])]
    index_array = da.stack(da.meshgrid(*indices,indexing='ij'), axis=-1).rechunk({-1:-1})

    @da.as_gufunc(signature=f"(i,j),({arr.ndim-2})->({arr.ndim-2})", output_dtypes=int, vectorize=True)
    def saveimg(image, index):
        imsave.save(fname, image)
        return index
    
    res = saveimg(arr,index_array)
    if compute == True:
        res.compute()
    else:
        return res


def run_cellpose_inference(
    im,
    model,
    diameter: float = 20,
    channels: list = [[0, 0]],
    flow_threshold: float = 0.5,
    cellprob_threshold: float = 0,
):
    """
    Perform cell segmentation using the Cellpose model.

    Args:
        im (numpy.ndarray): Input image to perform cell segmentation on.
        model: Cellpose model for segmentation.
        diameter (float, optional): Diameter parameter for cell segmentation. Defaults to 20.
        channels (list, optional): List of channels to use for segmentation. Defaults to [[0, 0]].
        flow_threshold (float, optional): Flow threshold for segmentation. Defaults to 0.5.
        cellprob_threshold (float, optional): Cell probability threshold. Defaults to 0.

    Returns:
        numpy.ndarray: Segmentation mask.

    Notes:
        - This function performs cell segmentation using the Cellpose model.
    """
    mask, _, _ = model.eval(
        [np.array(im)],
        diameter=diameter,
        channels=channels,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
    )
    return mask


def _segment_chunk(block, block_id, num_blocks, shift, **kwargs):
    # from squidpy
    if len(num_blocks) == 2:
        block_num = block_id[0] * num_blocks[1] + block_id[1]
    elif len(num_blocks) == 3:
        block_num = block_id[0] * (num_blocks[1] * num_blocks[2]) + block_id[1] * num_blocks[2]
    elif len(num_blocks) == 4:
        if num_blocks[-1] != 1:
            raise ValueError(
                f"Expected the number of blocks in the Z-dimension to be `1`, found `{num_blocks[-1]}`."
            )
        block_num = block_id[0] * (num_blocks[1] * num_blocks[2]) + block_id[1] * num_blocks[2]
    else:
        raise ValueError(f"Expected either `2`, `3` or `4` dimensional chunks, found `{len(num_blocks)}`.")

    labels = np.array(run_cellpose_inference(block, **kwargs)[0])
    mask = labels > 0
    labels[mask] = (labels[mask] << shift) | block_num

    return labels


def expand_labels(label_image, distance=1):
    from scipy.ndimage import distance_transform_edt
    def expand_labels_block(block):
        distances, nearest_label_coords = distance_transform_edt(
            block == 0, return_indices=True
        )
        labels_out = np.zeros_like(block)
        dilate_mask = distances <= distance
        masked_nearest_label_coords = [
            dimension_indices[dilate_mask]
            for dimension_indices in nearest_label_coords
        ]
        nearest_labels = block[tuple(masked_nearest_label_coords)]
        labels_out[dilate_mask] = nearest_labels
        return labels_out

    if isinstance(label_image, da.Array):
        return da.map_blocks(expand_labels_block, label_image, dtype=np.uint64)
    else:
        return expand_labels_block(label_image)

def _cellpose_segment(im, args):
    """
    Wrapper for whole or tiled segmentation with cellpose.

    Args:
        args: Argument object containing input and output file paths and parameters.

    Raises:
        FileNotFoundError: If input or output directories do not exist.
    """

    try:
        from cellpose import models
    except ImportError:
        raise ImportError("'cellpose' could not be found. Please install with 'pip install cellpose'")

    if args.metadata != "" and not check_directory_exists(args.metadata):
        raise FileNotFoundError("Parent directory for the metadata does not exist")
    
    _num_workers = 1
    if args.num_workers > 0:
        _num_workers = args.num_workers

    # Load cellpose model (path or pretrained)
    # TODO: check GPU available in torch
    _gpu = False
    if args.device == 'cuda':
        _gpu = True

    if args.model in OPENST_MODEL_NAMES:
        _model_path = cache_model_path(args.model)
        model = models.CellposeModel(gpu=_gpu, pretrained_model=_model_path)
    elif args.model in models.MODEL_NAMES:
        model = models.Cellpose(gpu=_gpu, model_type=args.model).cp
    elif check_file_exists(args.model):
        model = models.CellposeModel(gpu=_gpu, pretrained_model=args.model)
    else:
        raise ValueError(f"Cellpose model {args.model} was not found")

    if args.chunked:
        logging.info("Loading images into chunks")
        # TODO: implement checking of input file (dimensions)
        im = im.rechunk({0: args.chunk_size, 1: args.chunk_size})
        shift = int(np.prod(im.numblocks) - 1).bit_length()
        
        # This dask chunked option is useful for limited GPU memory
        # we need to specify processes scheduler such that CUDA
        # does not throw memory access errors
        # TODO: implement expand_labels (mask dilation) with chunked processing
        def _mask_tissue_wrapper(im):
            return mask_tissue(im, 
                        keep_black_background=args.keep_black_background, 
                        mask_gaussian_blur=args.tissue_masking_gaussian_sigma,
                        return_hsv=False)
            
        with ProgressBar():
            with dask.config.set(scheduler='single-threaded', num_workers=_num_workers):
                if args.mask_tissue:
                    logging.info("Masking whole tissue from background")
                    im.map_blocks(_mask_tissue_wrapper, meta=np.array((), dtype=np.int32))
                
                logging.info("Segmenting & relabeling by chunks")
                mask_chunked = da.map_overlap(
                    _segment_chunk,
                    im,
                    dtype=np.uint64,
                    num_blocks=im.numblocks,
                    shift=shift,
                    drop_axis=im.ndim - 1,
                    **{"model": model, "diameter": args.diameter, "flow_threshold": args.flow_threshold, "cellprob_threshold": args.cellprob_threshold}
                )
                label_groups = label_adjacency_graph(mask_chunked, None, mask_chunked.max())
                new_labeling = connected_components_delayed(label_groups)
                mask_complete = relabel_blocks(mask_chunked, new_labeling)
                if args.dilate_px > 0:
                    mask_complete = expand_labels(mask_complete, distance=args.dilate_px)
                
                if args.outline_px > 0:
                    mask_complete = find_boundaries(mask_complete, connectivity=1, mode='inner', background=0)
                    mask_complete = da.where(mask_complete, mask_complete, 0)
    else:
        if args.mask_tissue:
            logging.info("Masking whole tissue from background")
            im = mask_tissue(im, 
                            keep_black_background=args.keep_black_background, 
                            mask_gaussian_blur=args.tissue_masking_gaussian_sigma,
                            return_hsv=False)

        logging.info("Segmenting & relabeling whole image")
        mask_complete = run_cellpose_inference(
            im,
            model,
            diameter=args.diameter,
            flow_threshold=args.flow_threshold,
            cellprob_threshold=args.cellprob_threshold,
        )[0]

        if args.dilate_px > 0:
            mask_complete = expand_labels(mask_complete, distance=args.dilate_px)

        if args.outline_px > 0:
            mask_complete_bdy = find_boundaries(mask_complete, connectivity=1, mode='inner', background=0)
            mask_complete_bdy = expand_labels(mask_complete_bdy, distance=args.outline_px)
            mask_complete = np.where(mask_complete_bdy, mask_complete, 0)

        mask_complete = measure.label(mask_complete)

    return im, mask_complete

def _load_image(args):
    if args.h5_in != '':
        check_file_exists(args.h5_in)
        adata = h5py.File(args.h5_in, 'r+')

        if args.rna_segment:
            im, _ = create_unpaired_pseudoimage(adata, 
                                      args.rna_segment_spatial_coord_key,
                                      args.rna_segment_input_resolution,
                                      args.rna_segment_render_scale,
                                      args.rna_segment_render_sigma,
                                      args.rna_segment_output_resolution)
            logging.info("Will perform segmentation on a pseudoimage")
        else:
            im = adata[args.image_in]
            if args.chunked:
                im = da.from_array(im)
            else:
                im = im[:]
            logging.info("Will perform segmentation on a staining image")

    else:
        check_file_exists(args.image_in)
        if not check_directory_exists(args.mask_out):
            raise FileNotFoundError("Parent directory for --mask-out does not exist")
        if args.chunked:
            im = dask_image.imread.imread(args.image_in)[0] # to have 3 dimensions
        else:
            im = np.array(Image.open(args.image_in))
        
        adata = None
        
    return adata, im

def _save_mask(adata, im, mask_complete, args):
    # Transpose the image, so the axes are cxy
    if args.chunked:
        if args.h5_in != "":
            if args.mask_out in adata:
                logging.warn(f"The object {args.mask_out} will be removed from the h5py file")
                del adata[args.mask_out]

            dset = adata.create_dataset(args.mask_out, shape=mask_complete.shape,
                                            dtype=mask_complete.dtype)  
            logging.info(f'Saving mask to adata in {args.mask_out}')
            da.store(mask_complete, dset)

            logging.info(f'Relabeling mask in {args.mask_out}')
            adata[args.mask_out][...] = measure.label(adata[args.mask_out])
            
        else:
            # Save large image mask as zarr
            store = parse_url(args.mask_out, mode="w").store
            root = zarr.group(store=store)
            labels_grp = root.create_group('labels')
            logging.info(f'Saving mask to separate file in {args.mask_out}')
            write_image(im.transpose(2, 0, 1), group=root, compute=True, axes=['c', 'x', 'y'])
            write_image(mask_complete, group=labels_grp, compute=True, axes=['x', 'y'])
    else:
        dtype = np.uint8
        if mask_complete.max() >= (2**16 - 1):
            dtype = np.uint16
        elif mask_complete.max() >= (2**32 - 1):
            dtype = np.uint32
        elif mask_complete.max() >= (2**64 - 1):
            dtype = np.uint64

        if args.h5_in != "":
            if args.mask_out in adata:
                logging.warn(f"The object {args.mask_out} will be removed from the h5py file")
                del adata[args.mask_out]

            logging.info(f'Saving mask to adata in {args.mask_out}')
            adata[args.mask_out] = mask_complete
        else:
            logging.info(f'Saving mask to separate file in {args.mask_out}')
            Image.fromarray(mask_complete.astype(dtype)).save(args.mask_out)

def _run_segment(args):
    # Set maximum image size to support large HE (if not chunked)
    Image.MAX_IMAGE_PIXELS = args.max_image_pixels

    adata, im = _load_image(args)
    im, mask_complete = _cellpose_segment(im, args)
    _save_mask(adata, im, mask_complete, args)

if __name__ == "__main__":
    from openst.cli import get_segment_parser
    args = get_segment_parser().parse_args()
    _run_segment()
