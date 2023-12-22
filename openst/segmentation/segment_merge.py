import argparse
import h5py
import logging
import dask
import dask_image
import dask_image.imread
import dask.array as da
from dask.diagnostics import ProgressBar

import numpy as np
from PIL import Image
from skimage import measure
from ome_zarr.io import parse_url
from ome_zarr.writer import write_image
import zarr
from openst.utils.file import check_directory_exists, check_file_exists


def get_segment_merge_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="segmentation of Open-ST imaging data with cellpose",
        allow_abbrev=False,
        add_help=False,
    )
    parser.add_argument(
        "--mask-in",
        type=str,
        required=True,
        nargs=2,
        help="Path to the input segmentation masks - two of them!",
    )
    parser.add_argument(
        "--mask-out",
        type=str,
        required=True,
        help="Path (file or h5) where the merged mask will be saved",
    )
    parser.add_argument(
        "--adata",
        type=str,
        default='',
        help="When specified, masks are loaded from adata (--mask-in), and segmentation is saved there (to --mask-out)",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=512,
        help="When prediction of the mask runs in separate chunks, this is the chunk square size (in pixels)",
    )
    parser.add_argument(
        "--chunked",
        action="store_true",
        help="When specified, segmentation is computed at non-overlapping chunks of size '--chunk-size'",
    )
    parser.add_argument(
        "--max-image-pixels",
        type=int,
        default=933120000,
        help="Upper bound for number of pixels in the images (prevents exception when opening very large images)",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        help="Number of parallel workers when --chunked is specified",
        required=False,
        default=-1,
    )
    return parser


def setup_segment_merge_parser(parent_parser):
    """setup_segment_merge_parser"""
    parser = parent_parser.add_parser(
        "segment_merge",
        help="merge segmentations from Open-ST data",
        parents=[get_segment_merge_parser()],
    )
    parser.set_defaults(func=_run_segment_merge)

    return parser


def _segment_merge(mask_a, mask_b):
    mask_a = mask_a.astype(np.uint64)
    mask_b = mask_b.astype(np.uint64)
    mask_b = measure.label(mask_b)
    mask_a = measure.label(mask_a)
    mask_b[mask_a != 0] = 0
    mask_a = np.where(mask_a != 0, mask_a + mask_b.max(), 0)
    mask_out = mask_a + mask_b

    return mask_out

def _run_segment_merge(args):
    """
    Wrapper for segmentation merging

    Args:
        args: Argument object containing input and output file paths and parameters.

    Raises:
        FileNotFoundError: If input or output directories do not exist.
    """
    logging.info("openst segment merge; running with parameters:")
    logging.info(args.__dict__)

    # Check input and output data
    mask_a, mask_b = None, None
    if args.adata != '':
        check_file_exists(args.adata)
        adata = h5py.File(args.adata, 'r+')
        mask_a = adata[args.mask_in[0]]
        mask_b = adata[args.mask_in[1]]
        if args.chunked:
            mask_a = da.from_array(mask_a)
            mask_b = da.from_array(mask_b)
        else:
            mask_a = mask_a[:]
            mask_b = mask_b[:]
    else:
        check_file_exists(args.mask_in[0])
        check_file_exists(args.mask_in[1])
        if not check_directory_exists(args.mask_out):
            raise FileNotFoundError("Parent directory for --mask-out does not exist")
    
    _num_workers = 1
    if args.num_workers > 0:
        _num_workers = args.num_workers

    Image.MAX_IMAGE_PIXELS = args.max_image_pixels

    if args.chunked:
        logging.info("Loading images into chunks")
        if mask_a is None or mask_b is None:
            mask_a = dask_image.imread.imread(args.mask_in[0])[0]
            mask_b = dask_image.imread.imread(args.mask_in[1])[0]

        # rechunk both images
        mask_a = mask_a.rechunk({0: args.chunk_size, 1: args.chunk_size})
        mask_b = mask_b.rechunk({0: args.chunk_size, 1: args.chunk_size})

        logging.info("Merging segmentation masks by chunks")
        with ProgressBar():
            with dask.config.set(scheduler='single-threaded', num_workers=_num_workers):                
                _mask_out = _segment_merge(mask_a, mask_b)

                if args.adata:
                    if args.mask_out in adata:
                        logging.warn(f"The object {args.mask_out} will be removed from the h5py file")
                        del adata[args.mask_out]

                    dset = adata.create_dataset(args.mask_out, shape=_mask_out.shape,
                                                    dtype=_mask_out.dtype)  
                    logging.info(f'Saving mask to adata in {args.mask_out}')
                    da.store(_mask_out, dset)
                    
                else:
                    print(_mask_out)
                    store = parse_url(args.mask_out, mode="w").store
                    root = zarr.group(store=store)
                    labels_grp = root.create_group('labels')
                    logging.info(f'Saving mask to separate file in {args.mask_out}')
                    write_image(_mask_out.transpose(2, 0, 1), group=root, compute=True, axes=['c', 'x', 'y'])
                    write_image(_mask_out, group=labels_grp, compute=True, axes=['x', 'y'])
    else:
        if mask_a is None or mask_b is None:
            mask_a = np.array(Image.open(args.mask_in[0]))
            mask_b = np.array(Image.open(args.mask_in[1]))

        _mask_out = _segment_merge(mask_a, mask_b)

        if args.adata:
            if args.mask_out in adata:
                logging.warn(f"The object {args.mask_out} will be removed from the h5py file")
                del adata[args.mask_out]

            logging.info(f'Saving mask to adata in {args.mask_out}')
            adata[args.mask_out] = _mask_out
        else:
            logging.info(f'Saving mask to separate file in {args.mask_out}')
            Image.fromarray(_mask_out.astype(np.uint8)).save(args.mask_out)


if __name__ == "__main__":
    args = get_segment_merge_parser().parse_args()
    _run_segment_merge(args)
