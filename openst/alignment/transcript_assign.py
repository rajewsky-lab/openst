import logging

import numpy as np
import pandas as pd
from anndata import read_h5ad
from PIL import Image
from skimage import measure

from openst.utils.file import (check_directory_exists, check_file_exists,
                               load_properties_from_adata)
from openst.utils.spacemake import reassign_indices_adata


def assert_valid_mask(im):
    if len(im.shape) < 2:
        raise ValueError("The specified image mask is not two-dimensional")
    elif len(im.shape) > 2:
        if im.shape[:-1] != 1:
            raise ValueError("A 3D image should be XYC, where C=1 for a segmentation mask")


def transfer_segmentation(adata_transformed_coords, props_filter):
    joined_coordinates = np.array([props_filter["centroid-0"], props_filter["centroid-1"]]).T
    joined_coordinates = np.vstack([np.array([0, 0]), joined_coordinates])

    cell_ID_merged = np.array(props_filter["label"])
    cell_ID_merged = np.hstack([np.array([0]), cell_ID_merged])

    adata_by_cell = reassign_indices_adata(
        adata_transformed_coords,
        np.array(adata_transformed_coords.obs["cell_ID"]),
        joined_coordinates,
        cell_ID_merged,
    )

    _missing_uns_keys = set(list(adata_by_cell.uns.keys()) + list(adata_transformed_coords.uns.keys()))
    _missing_uns_keys = set(list(adata_transformed_coords.uns.keys())).intersection(_missing_uns_keys)

    for _missing_key in _missing_uns_keys:
        adata_by_cell.uns[_missing_key] = adata_transformed_coords.uns[_missing_key]

    return adata_by_cell


def subset_adata_to_mask(mask, adata, spatial_key: str = 'spatial'):
    # Subset adata to the valid coordinates from the mask
    adata = adata[(adata.obsm[spatial_key][:, 0] >= 0) & 
                  (adata.obsm[spatial_key][:, 1] >= 0) &
                  (adata.obsm[spatial_key][:, 0] < mask.shape[0]) & 
                  (adata.obsm[spatial_key][:, 1] < mask.shape[1])].copy()

    # Subset the labels to those in the mask
    labels = mask[adata.obsm[spatial_key][:, 0].astype(int), adata.obsm[spatial_key][:, 1].astype(int)]

    # Set first label as zero for background
    labels[0] = 0

    # Assign label as cell_ID
    adata.obs["cell_ID"] = labels

    # Get centroid and label ID from mask
    props = measure.regionprops_table(mask, properties=["label", "centroid"])
    props = pd.DataFrame(props)

    props_filter = props[props.label.isin(np.unique(adata.obs["cell_ID"]))]
    return adata, props_filter


def shuffle_umi(adata, spatial_key='spatial'):
    from scipy.sparse import csr_array, csc_matrix
    import anndata as ad

    X_0_repeat = np.repeat(adata.X.nonzero()[0], adata.X.data.astype(int))
    X_1_repeat = np.repeat(adata.X.nonzero()[1], adata.X.data.astype(int))
    tile_id = np.repeat(adata.obs['tile_id'].astype(str).values, np.array(adata.X.sum(axis=1)).flatten().astype(int))
    idx_range = np.arange(len(X_0_repeat))
    loc_random = np.random.randint(0, len(adata.obsm[spatial_key]), len(X_1_repeat))
    # obsm_spatial_expanded = tiles_transformed_coords_refined_concatenated[loc_random]
    obsm_spatial_expanded = adata.obsm[spatial_key][loc_random]
    X_repeat = csc_matrix(csr_array((np.ones_like(X_0_repeat), (idx_range, X_1_repeat))))
    adata_shuffled = ad.AnnData(X=X_repeat)
    adata_shuffled.obs_names = (pd.Series(idx_range.astype(str)) + ":" + pd.Series(tile_id)).values
    adata_shuffled.obsm[spatial_key] = obsm_spatial_expanded

    for col in adata.obs.columns:
        adata_shuffled.obs[col] = 1

    adata_shuffled.var_names = adata.var_names

    return adata_shuffled


def _run_transcript_assign(args):
    """_run_transcript_assign
    
    This one uses AnnData instead of h5py, so the obsm/ and uns/ keys must be parsed 
    """
    # TODO: load with dask if it is too large

    Image.MAX_IMAGE_PIXELS = args.max_image_pixels

    check_file_exists(args.h5_in)

    if args.mask_from_file:
        check_file_exists(args.mask_in)

    if not check_directory_exists(args.h5_out):
        raise FileNotFoundError("Parent directory for --h5-out does not exist")

    if args.metadata != "" and not check_directory_exists(args.metadata):
        raise FileNotFoundError("Parent directory for the metadata does not exist")
    
    if not args.spatial_key.startswith("obsm/"):
        raise ValueError("The '--spatial-key' must start with 'obsm/'")
    
    spatial_key = "/".join(args.spatial_key.split("/")[1:])

    if not args.mask_from_file:
        logging.info(f"Loading image mask from Open-ST h5 object at '{args.mask_in}'")
        mask = load_properties_from_adata(args.h5_in, [args.mask_in])[args.mask_in]
    else:
        logging.info(f"Loading image mask from file at '{args.mask_in}'")
        mask = np.array(Image.open(args.mask_in))

    assert_valid_mask(mask)

    logging.info("Parsing Open-ST h5 object as AnnData")
    adata = read_h5ad(args.h5_in)

    if args.shuffle_umi:
        logging.info("Spatially shuffling UMIs")
        adata = shuffle_umi(adata, spatial_key=spatial_key)

    logging.info("Subsetting Open-ST AnnData coordinates to mask")
    adata, props_filter = subset_adata_to_mask(mask, adata, spatial_key)

    logging.info("Assigning transcripts to segmented cells")
    adata_by_cell = transfer_segmentation(adata, props_filter)

    logging.info(f"Writing Open-ST AnnData by segmented cells to {args.h5_out}")
    adata_by_cell.write_h5ad(args.h5_out)


if __name__ == "__main__":
    from openst.cli import get_transcript_assign_parser
    args = get_transcript_assign_parser().parse_args()
    _run_transcript_assign(args)
