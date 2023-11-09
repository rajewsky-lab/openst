import argparse
import logging

import numpy as np
import pandas as pd
from anndata import read_h5ad
from PIL import Image
from skimage import measure

from openst.utils.file import (check_directory_exists, check_file_exists,
                               load_properties_from_adata)
from openst.utils.spacemake import reassign_indices_adata


def get_transcript_assign_parser():
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False,
        description="openst transfer of transcripts to single cells using a pairwise-aligned segmentation mask",
    )

    parser.add_argument(
        "--adata",
        type=str,
        help="path to previously aligned spatial.h5ad AnnData file",
        required=True,
    )

    parser.add_argument(
        "--mask-in-adata",
        default=False,
        action="store_true",
        help="When specified, the image mask is loaded from the adata, at the internal path specified by '--mask'",
    )

    parser.add_argument(
        "--shuffle-umi",
        default=False,
        action="store_true",
        help="When specified, UMI locations will be shuffled. This can be used as a baseline for feature selection.",
    )

    parser.add_argument(
        "--mask",
        type=str,
        help="""path to image mask; must be in same coordinates as the obsm['spatial'] in the AnnData file.
        If --mask-in-adata, it is a path within the h5ad file""",
        required=True,
    )

    parser.add_argument(
        "--spatial-key",
        type=str,
        help="""obsm dataset for the aligned coordinates, e.g. 'spatial_pairwise_aligned_coarse'  or
        'spatial_pairwise_aligned_fine' if the data has been aligned with openst pairwise_aligner""",
        required=True,
    )

    parser.add_argument(
        "--output",
        type=str,
        help="path and filename for output file that will be generated",
        required=True,
    )

    parser.add_argument(
        "--max-image-pixels",
        type=int,
        default=933120000,
        help="Upper bound for number of pixels in the images (prevents exception when opening very large images)",
    )
    parser.add_argument(
        "--metadata-out",
        type=str,
        default="",
        help="""Path where the metadata will be stored.
        If not specified, metadata is not saved.
        Warning: a report (via openst report) cannot be generated without metadata!""",
    )

    return parser


def setup_transcript_assign_parser(parent_parser):
    """setup_transcript_assign_parser"""
    parser = parent_parser.add_parser(
        "transcript_assign",
        help="assign transcripts into previously aligned segmentation mask",
        parents=[get_transcript_assign_parser()],
    )
    parser.set_defaults(func=_run_transcript_assign)

    return parser


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
    adata = adata[(adata.obsm[spatial_key][:, 0] <= mask.shape[0]) & (adata.obsm[spatial_key][:, 1] <= mask.shape[1])].copy()

    # Subset the labels to those in the mask
    labels = mask[adata.obsm[spatial_key][:, 0].astype(int), adata.obsm[spatial_key][:, 1].astype(int)]

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
    """_run_transcript_assign."""
    # TODO: load with dask if it is too large
    logging.info("openst spatial transcriptomics stitching; running with parameters:")
    logging.info(args.__dict__)

    Image.MAX_IMAGE_PIXELS = args.max_image_pixels

    check_file_exists(args.adata)

    if not args.mask_in_adata:
        check_file_exists(args.mask)

    if not check_directory_exists(args.output):
        raise FileNotFoundError("Parent directory for --output does not exist")

    if args.metadata_out != "" and not check_directory_exists(args.metadata_out):
        raise FileNotFoundError("Parent directory for the metadata does not exist")

    logging.info("Loading image mask")
    if args.mask_in_adata:
        mask = load_properties_from_adata(args.adata, [args.mask])[args.mask]
    else:
        mask = np.array(Image.open(args.mask))

    assert_valid_mask(mask)

    logging.info("Loading adata")
    adata = read_h5ad(args.adata)

    if args.shuffle_umi:
        logging.info("Spatially shuffling UMIs")
        adata = shuffle_umi(adata, spatial_key=args.spatial_key)

    logging.info("Subsetting adata coordinates to mask")
    adata, props_filter = subset_adata_to_mask(mask, adata, args.spatial_key)

    logging.info("Assigning transcripts to cells in mask")
    adata_by_cell = transfer_segmentation(adata, props_filter)

    logging.info(f"Writing output to {args.output}")
    adata_by_cell.write_h5ad(args.output)


if __name__ == "__main__":
    args = get_transcript_assign_parser().parse_args()
    _run_transcript_assign(args)
