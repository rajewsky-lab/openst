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


def transfer_segmentation(adata_transformed_coords, label_image, props_filter):
    joined_coordinates = np.array([props_filter["centroid-0"], props_filter["centroid-1"]]).T
    joined_coordinates = np.vstack([np.array([0, 0]), joined_coordinates])

    cell_ID_merged = np.array(props_filter["label"])
    cell_ID_merged = np.hstack([np.array([0]), cell_ID_merged])

    adata_by_cell = reassign_indices_adata(
        adata_transformed_coords,
        np.array(adata_transformed_coords.obs["cell_ID"]),
        joined_coordinates,
        label_image,
        cell_ID_merged,
    )

    spatial_units_obs_names_dict = {}

    for sn in adata_by_cell.uns["spatial_units_obs_names"]:
        bc, tile = sn.split(":")
        if tile in spatial_units_obs_names_dict.keys():
            spatial_units_obs_names_dict[tile] += [bc]
        else:
            spatial_units_obs_names_dict[tile] = [bc]

    for sn in adata_by_cell.uns["spatial_units_obs_names"]:
        bc, tile = sn.split(":")
        if tile in spatial_units_obs_names_dict.keys():
            spatial_units_obs_names_dict[tile] += [bc]
        else:
            spatial_units_obs_names_dict[tile] = [bc]

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

    # we transpose the mask to have the same coordinates
    mask = mask.T

    logging.info("Loading adata")
    adata = read_h5ad(args.adata)

    logging.info("Subsetting adata coordinates to mask")
    adata, props_filter = subset_adata_to_mask(mask, adata, args.spatial_key)

    logging.info("Assigning transcripts to cells in mask")
    adata_by_cell = transfer_segmentation(adata, mask, props_filter)

    logging.info(f"Writing output to {args.output}")
    adata_by_cell.write_h5ad(args.output)


if __name__ == "__main__":
    args = get_transcript_assign_parser().parse_args()
    _run_transcript_assign(args)
