"""
Automatic Pairwise Alignment of Spatial Transcriptomics and Imaging Data (open-ST)

Author: Daniel León-Periñán @ N.Rajewsky Lab (BIMSB)
Date: August 30, 2023
Version: 0.1.0

Description:
This script performs automatic pairwise alignment between spatial transcriptomics data and imaging data.
It aligns the spatial transcriptomics data onto the imaging data to enable spatial comparison and analysis.

Usage example:
openst pairwise_aligner --image-in input_image.jpg --h5-in input_data.h5ad --h5-out aligned_data.h5ad
                     --metadata-out metadata.pkl --save-image-in-h5 --rescale-factor-coarse 20
                     --rescale-factor-fine 5 --tissue-masking-gaussian-sigma 5 --fine-registration-gaussian-sigma 2
                     --keep-black-background --threshold-counts-coarse 1 --threshold-counts-fine 0
                     --pseudoimage-size-coarse 4000 --pseudoimage-size-fine 6000 --ransac-coarse-min-samples 3
                     --ransac-coarse-residual-threshold 2 --ransac-coarse-max-trials 10000 --ransac-fine-min-samples 10
                     --ransac-fine-residual-threshold 2 --ransac-fine-max-trials 10000 --max-image-pixels 933120000
                     --feature-matcher 'LoFTR' --n-threads 2 --fine-min-matches 50 --fiducial-model model.pt
"""

import argparse
import logging

import numpy as np
from anndata import read_h5ad
from PIL import Image
from skimage.transform import estimate_transform
from threadpoolctl import threadpool_limits

from openst.alignment.transformation import apply_transform
from openst.metadata.classes.pairwise_alignment import (
    AlignmentResult, PairwiseAlignmentMetadata)
from openst.utils.file import (check_adata_structure, check_directory_exists,
                               check_file_exists, load_properties_from_adata)
from openst.utils.pseudoimage import create_pseudoimage


def get_manual_pairwise_aligner_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="openst pairwise alignment of two-dimensional spatial transcriptomics and imaging data",
        allow_abbrev=False,
        add_help=False,
    )

    parser.add_argument(
        "--image-in",
        type=str,
        required=True,
        help="Path to the input image. This is treated as the 'destination' image during pairwise alignment",
    )
    parser.add_argument(
        "--h5-in",
        type=str,
        required=True,
        help="Path to the input h5ad file containing spatial coordinates",
    )
    parser.add_argument(
        "--metadata-out",
        type=str,
        default="",
        help="Path where the metadata will be stored. If not specified, metadata is not saved.",
    )
    parser.add_argument(
        "--h5-out",
        type=str,
        required=True,
        help="Path where the h5ad file will be saved after alignment",
    )
    parser.add_argument(
        "--save-image-in-h5",
        action="store_true",
        help="Whether the input image will be saved into the output h5ad file",
    )
    parser.add_argument(
        "--coarse",
        action="store_true",
        help="If selected, manual coarse alignment will be performed",
    )
    parser.add_argument(
        "--fine",
        action="store_true",
        help="If selected, fine (per tile) coarse alignment will be performed",
    )
    parser.add_argument(
        "--rescale-factor-coarse",
        type=int,
        default=20,
        help="Rescaling factor for the input image (1:factor), used during coarse pairwise alignment",
    )
    parser.add_argument(
        "--rescale-factor-fine",
        type=int,
        default=5,
        help="Rescaling factor for the input image (1:factor), used during fine pairwise alignment",
    )
    parser.add_argument(
        "--threshold-counts-coarse",
        type=int,
        default=1,
        help="""Only spatial coordinates with counts larger than this number
        will be kept for pseudoimage rendering during coarse alignment""",
    )
    parser.add_argument(
        "--threshold-counts-fine",
        type=int,
        default=0,
        help="""Only spatial coordinates with counts larger than this number
        will be kept for pseudoimage rendering during fine alignment""",
    )
    parser.add_argument(
        "--pseudoimage-size-coarse",
        type=int,
        default=4000,
        help="Size (in pixels) of the pseudoimage during coarse alignment.",
    )
    parser.add_argument(
        "--pseudoimage-size-fine",
        type=int,
        default=6000,
        help="Size (in pixels) of the pseudoimage during fine alignment.",
    )
    parser.add_argument(
        "--max-image-pixels",
        type=int,
        default=933120000,
        help="Upper bound for number of pixels in the images (prevents exception when opening very large images)",
    )
    return parser


def setup_manual_pairwise_aligner_parser(parent_parser):
    """setup_manual_pairwise_aligner_parser"""
    parser = parent_parser.add_parser(
        "manual_pairwise_aligner",
        help="openst manual pairwise alignment of spatial transcriptomics and imaging data",
        parents=[get_manual_pairwise_aligner_parser()],
    )
    parser.set_defaults(func=_run_manual_pairwise_aligner)

    return parser


def run_manual_registration(
    in_coords: np.ndarray,
    total_counts: np.ndarray,
    puck_id: np.ndarray,
    staining_image: np.ndarray,
    args,
) -> (np.ndarray, np.ndarray, PairwiseAlignmentMetadata):
    """
    Perform manual registration of spatial transcriptomics (STS) data with a staining image.

    Args:
        in_coords (np.ndarray): Input STS coordinates.
        total_counts (np.ndarray): Total UMI counts for each STS coordinate.
        puck_id: Identifier for each STS coordinate. During the fine registration,
                 this 'puck_id' is used to aggregate the coordinates into buckets that
                 are aligned separately. Recommended for flow-cell based STS.
        staining_image (np.ndarray): Staining image for registration.
        args: Namespace containing various registration parameters.

    Returns:
        tuple: A tuple containing four elements:
            - out_coords_output_coarse (np.ndarray): Registered STS coordinates after coarse registration
            - out_coords_output_fine (np.ndarray): Registered STS coordinates after fine registration
            - registered_staining_image (np.ndarray): Staining image after registration.
            - metadata (PairwiseAlignmentMetadata)
    """
    # Create output objects
    out_coords_output_fine = np.zeros_like(in_coords)
    metadata = PairwiseAlignmentMetadata(args)

    logging.info(f"Coarse registration, {len(in_coords)} coordinates")

    # Preparing images and preprocessing routines
    sts_coords = in_coords[total_counts > args.threshold_counts_coarse]
    staining_image_rescaled = staining_image[:: args.rescale_factor_coarse, :: args.rescale_factor_coarse]
    sts_pseudoimage = create_pseudoimage(sts_coords, args.pseudoimage_size_coarse, staining_image_rescaled.shape)

    # Feature matching
    # TODO: here, get the points
    in_mkpts0 = []
    in_mkpts1 = []
    sts_pseudoimage["pseudoimage"]
    staining_image_rescaled

    # Estimate and apply transform
    tform_points = estimate_transform("similarity", in_mkpts0, in_mkpts1)
    sts_coords_to_transform = sts_pseudoimage["coords_rescaled"] * sts_pseudoimage["rescaling_factor"]
    sts_coords_transformed = apply_transform(sts_coords_to_transform, tform_points, check_bounds=False)
    sts_coords_transformed = sts_coords_transformed * args.rescale_factor_coarse

    # Remove coordinates outside the image (for pseudoimage generation)
    _i_sts_coords_coarse_within_image_bounds = np.where(
        (sts_coords_transformed[:, 0] > 0)
        & (sts_coords_transformed[:, 0] < staining_image.shape[1])
        & (sts_coords_transformed[:, 1] > 0)
        & (sts_coords_transformed[:, 1] < staining_image.shape[0])
    )
    sts_coords_transformed = sts_coords_transformed[_i_sts_coords_coarse_within_image_bounds][:, :-1]

    # Apply transform to all coordinates & retransform back
    sts_coords_coarse = in_coords.copy()
    sts_coords_coarse -= sts_pseudoimage["offset_factor"]
    sts_coords_coarse = (
        (sts_coords_coarse / sts_pseudoimage["rescale_factor"]) * sts_pseudoimage["scale"]
    ) * sts_pseudoimage["rescaling_factor"]

    sts_coords_coarse = apply_transform(sts_coords_coarse, tform_points, check_bounds=False)
    sts_coords_coarse = sts_coords_coarse * args.rescale_factor_coarse
    sts_coords_coarse = sts_coords_coarse[:, :-1]
    out_coords_output_coarse = sts_coords_coarse.copy()

    # Saving alignment results here
    _align_result = AlignmentResult(
        im_0=None,
        im_1=sts_pseudoimage["pseudoimage"],
        transformation_matrix=tform_points.params,
        ransac_results=None,
        sift_results=None,
        keypoints0=in_mkpts0,
        keypoints1=in_mkpts1,
    )
    metadata.add_alignment_result(_align_result)

    # Finish here if only coarse registration was selected to run
    if args.only_coarse:
        return (
            out_coords_output_coarse,
            None,
            staining_image,
            metadata,
        )

    # STAGE 2: fine registration per tile
    logging.info("Running manual fine registration (per tile)")

    # Collect tile identifiers
    tile_codes = np.unique(puck_id.codes)

    # Apply scaling to input image again, for fine registration
    staining_image_rescaled = staining_image[:: args.rescale_factor_fine, :: args.rescale_factor_fine]

    for tile_code in tile_codes:
        # Create a pseudoimage
        _t_puck_id = puck_id.codes == tile_code
        _t_valid_coords = (
            puck_id[(total_counts > args.threshold_counts_coarse)].codes[_i_sts_coords_coarse_within_image_bounds]
            == tile_code
        )

        logging.info(f"Registering tile {tile_code} with {len(_t_valid_coords)} coordinates")

        _t_sts_pseudoimage = create_pseudoimage(
            sts_coords_transformed[:, ::-1],  # we need to flip these coordinates
            args.pseudoimage_size_fine,
            staining_image_rescaled.shape,
            _t_valid_coords,
            recenter=False,
            rescale=True,
            values=None,
        )

        # Axis limits to crop both modalities to tile region
        _t_sts_coords_to_transform = _t_sts_pseudoimage["coords_rescaled"] * _t_sts_pseudoimage["rescaling_factor"]

        min_lim, max_lim = _t_sts_coords_to_transform[_t_valid_coords].min(axis=0).astype(
            int
        ), _t_sts_coords_to_transform[_t_valid_coords].max(axis=0).astype(int)
        x_min, y_min = min_lim
        x_max, y_max = max_lim

        # Finding matches between modalities
        # TODO: get the manual coordinates here
        _t_mkpts0, _t_mkpts1 = [], []
        _t_sts_pseudoimage["pseudoimage"][x_min:x_max, y_min:y_max]
        staining_image_rescaled[x_min:x_max, y_min:y_max]

        # Compute similarity matrix and compute point transformation
        _t_tform_points = estimate_transform("similarity", _t_mkpts0, _t_mkpts1)

        # Apply the same transformation to the tiles
        _t_sts_coords_fine_to_transform = sts_coords_coarse[_t_puck_id] / args.rescale_factor_fine
        _t_sts_coords_fine_to_transform = (_t_sts_coords_fine_to_transform - np.array([[y_min, x_min]]))[:, ::-1]

        _t_sts_coords_fine_transformed = apply_transform(
            _t_sts_coords_fine_to_transform, _t_tform_points, check_bounds=True
        )

        # Rescale points to original HE dimensions
        _t_sts_coords_fine_transformed = _t_sts_coords_fine_transformed[:, :-1] + np.array([[y_min, x_min]])
        _t_sts_coords_fine_transformed = _t_sts_coords_fine_transformed * args.rescale_factor_fine

        out_coords_output_fine[_t_puck_id] = _t_sts_coords_fine_transformed

        # Saving alignment results here (only when passed)
        _align_result = AlignmentResult(
            im_0=None,
            im_1=sts_pseudoimage["pseudoimage"],
            transformation_matrix=tform_points.params,
            ransac_results=None,
            sift_results=None,
            keypoints0=in_mkpts0,
            keypoints1=in_mkpts1,
        )
        metadata.add_alignment_result(_align_result)

    return (
        out_coords_output_coarse,
        out_coords_output_fine,
        staining_image,
        metadata,
    )


def _run_manual_pairwise_aligner(args):
    logging.info("open-ST pairwise alignment; running with parameters:")
    logging.info(args.__dict__)

    # Check input and output data
    check_file_exists(args.h5_in)
    check_adata_structure(args.h5_in)
    check_file_exists(args.image_in)

    if not check_directory_exists(args.h5_out):
        raise FileNotFoundError("Parent directory for --h5-out does not exist")

    if args.metadata_out != "" and not check_directory_exists(args.metadata_out):
        raise FileNotFoundError("Parent directory for the metadata does not exist")

    if args.h5_out == args.h5_in:
        raise ValueError("The path to the output file cannot be the same as the input file")

    # Set maximum image size to support large HE
    Image.MAX_IMAGE_PIXELS = args.max_image_pixels

    # Loading the spatial transcriptomics data
    sts = load_properties_from_adata(args.h5_in, properties=["obsm/spatial", "obs/total_counts", "obs/puck_id"])

    # Loading image data
    staining_image = np.array(Image.open(args.image_in))

    # Running registration
    sts_aligned_coarse, sts_aligned_fine, staining_image_aligned, metadata = run_manual_registration(
        sts["obsm/spatial"],
        sts["obs/total_counts"],
        sts["obs/puck_id"],
        staining_image,
        args,
    )

    # Saving the metadata (for QC)
    if args.metadata_out != "":
        metadata.render_images()
        metadata.save_json(args.metadata_out)

    # Exporting the data
    logging.info(f"Loading {args.h5_in}")
    adata = read_h5ad(args.h5_in)
    adata.obsm["spatial_pairwise_aligned_coarse"] = sts_aligned_coarse
    if sts_aligned_fine is not None:
        adata.obsm["spatial_pairwise_aligned_fine"] = sts_aligned_fine

    if args.save_image_in_h5:
        adata.uns["spatial_pairwise_aligned"] = {}
        adata.uns["spatial_pairwise_aligned"]["staining_image"] = staining_image
        adata.uns["spatial_pairwise_aligned"]["staining_image_transformed"] = staining_image_aligned

    logging.info(f"Writing output to {args.h5_out}")
    adata.write(args.h5_out)
    logging.info(f"Output {args.h5_out} file was written")


if __name__ == "__main__":
    args = get_manual_pairwise_aligner_parser().parse_args()
    with threadpool_limits(limits=args.n_threads):
        _run_manual_pairwise_aligner(args)
