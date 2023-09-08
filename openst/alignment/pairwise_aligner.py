#!/usr/bin/env python3

"""
Automatic Pairwise Alignment of Spatial Transcriptomics and Imaging Data (open-ST)

Author: Daniel León-Periñán @ N.Rajewsky Lab (BIMSB)
Date: August 30, 2023
Version: 0.1.0

Description:
This script performs automatic pairwise alignment between spatial transcriptomics data and imaging data.
It aligns the spatial transcriptomics data onto the imaging data to enable spatial comparison and analysis.

Usage example:
python pairwise_aligner.py --image-in input_image.jpg --h5-in input_data.h5ad --h5-out aligned_data.h5ad
                     --metadata-out metadata.pkl --save-image-in-h5 --rescale-factor-coarse 20
                     --rescale-factor-fine 5 --tissue-masking-gaussian-sigma 5 --fine-registration-gaussian-sigma 2
                     --keep-black-background --threshold-counts-coarse 1 --threshold-counts-fine 0
                     --pseudoimage-size-coarse 4000 --pseudoimage-size-fine 6000 --ransac-coarse-min-samples 3
                     --ransac-coarse-residual-threshold 2 --ransac-coarse-max-trials 10000 --ransac-fine-min-samples 10
                     --ransac-fine-residual-threshold 2 --ransac-fine-max-trials 10000 --max-image-pixels 933120000
                     --feature-matcher 'SIFT' --n-threads 2
"""

import argparse
import logging
import numpy as np

from anndata import read_h5ad
from itertools import product
from PIL import Image
from scipy import ndimage
from skimage.color import rgb2hsv
from skimage.exposure import equalize_adapthist
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import ransac
from skimage.transform import (SimilarityTransform, estimate_transform, rotate)
from threadpoolctl import threadpool_limits

from openst.file import save_pickle, check_file_exists, check_adata_structure, load_properties_from_adata, check_directory_exists
from openst.alignment import feature_matching
from openst.alignment.metadata import PairwiseAlignmentMetadata, AlignmentResult
from openst.utils.image import create_pseudoimage

ROTATION_INVARIANT = ['SIFT']

def get_pairwise_aligner_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="open-ST pairwise alignment of two-dimensional spatial transcriptomics and imaging data.",
        allow_abbrev=False,
        add_help=False
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
        "--mask-tissue",
        action="store_true",
        help="Tissue (imaging modality) is masked from the background for the feature detection",
    )
    parser.add_argument(
        "--only-coarse",
        action="store_true",
        help="If selected, only the coarse alignment stage will run",
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
        "--tissue-masking-gaussian-sigma",
        type=int,
        default=5,
        help="The gaussian blur sigma used during the isolation of the tissue on the HE (preprocessing)",
    )
    parser.add_argument(
        "--fine-registration-gaussian-sigma",
        type=int,
        default=2,
        help="Gaussian blur used on all modalities during fine registration",
    )
    parser.add_argument(
        "--keep-black-background",
        action="store_true",
        help="Whether to set the background of the imaging modalities to white, after tissue masking",
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
        "--ransac-coarse-min-samples",
        type=int,
        default=3,
        help="'min_samples' parameter of RANSAC, during coarse registration",
    )
    parser.add_argument(
        "--ransac-coarse-residual-threshold",
        type=float,
        default=2,
        help="'residual_threshold' parameter of RANSAC, during coarse registration",
    )
    parser.add_argument(
        "--ransac-coarse-max-trials",
        type=int,
        default=10000,
        help="'max_trials' parameter of RANSAC, during coarse registration",
    )
    parser.add_argument(
        "--ransac-fine-min-samples",
        type=int,
        default=10,
        help="'min_samples' parameter of RANSAC, during fine registration",
    )
    parser.add_argument(
        "--ransac-fine-residual-threshold",
        type=float,
        default=2,
        help="'residual_threshold' parameter of RANSAC, during fine registration",
    )
    parser.add_argument(
        "--ransac-fine-max-trials",
        type=int,
        default=10000,
        help="'max_trials' parameter of RANSAC, during fine registration",
    )
    parser.add_argument(
        "--max-image-pixels",
        type=int,
        default=933120000,
        help="'max_trials' parameter of RANSAC, during fine registration",
    )
    parser.add_argument(
        "--n-threads",
        type=int,
        default=1,
        help="Number of CPU threads for parallel processing",
    )
    parser.add_argument(
        "--feature-matcher",
        type=str,
        default="LoFTR",
        choices=["LoFTR", "SIFT"],
        help="Feature matching algorithm",
    )
    parser.add_argument(
        "--fine-min-matches",
        type=int,
        default=50,
        help="Minimum number of matching keypoints between modalities during fine alignment",
    )
    return parser


def setup_pairwise_aligner_parser(parent_parser):
    """setup_pairwise_aligner_parser"""
    parser = parent_parser.add_parser(
        "pairwise_aligner",
        help="simulate pairwise_aligner on PASTA-seq data",
        parents=[get_pairwise_aligner_parser()],
    )
    parser.set_defaults(func=_run_pairwise_aligner)

    return parser


def prepare_image_for_feature_matching(
    image: np.ndarray,
    gaussian_blur: float = 0,
    flip: list = [1, 1],
    rotation: float = 0,
    crop: list = [9, None, 0, None],
    mask_tissue: bool = False,
    keep_black_background: bool = False,
    mask_gaussian_blur: float = 5,
):
    """
    Prepare an image for feature matching by applying a series of image processing steps.

    Args:
        image (np.ndarray): Input RGB image to be prepared.
        gaussian_blur (float, optional): Standard deviation for Gaussian blurring applied after processing.
        flip (list, optional): List indicating whether to flip the image along the x and y axes.
        rotation (float, optional): float indicating the rotation angle in degrees in counter-clockwise direction.
        crop (list, optional): List specifying cropping range [x_min, x_max, y_min, y_max] for the prepared image.
        mask_tissue (bool, optional): If True, background is removed with saturation -> blur -> otsu
        keep_black_background (bool, optional): If True, keeps background pixels as black. Only if mask_tissue=True
        mask_gaussian_blur (float, optional): Standard deviation for Gaussian blurring of the saturation channel.

    Returns:
        list: A list containing prepared images after applying the specified processing steps.
            - Each element of the list corresponds to a processed channel of the input image.

    Notes:
        - This function prepares the input image for feature matching by applying a series of processing steps.
        - The steps include thresholding, binary fill, background removal, equalization, blurring, and flipping.
        - The 'mask_gaussian_blur' parameter controls the amount of blurring applied to the saturation channel.
        - The 'keep_black_background' parameter toggles background pixel replacement with white.
        - The 'gaussian_blur' parameter controls the amount of blurring applied after processing.
        - The 'flip' parameter specifies flipping along the x and y axes using a list of factors [x_flip, y_flip].
        - The 'crop' parameter defines the cropping range for the prepared image.

    Raises:
        ValueError: If provided arguments have invalid values.
        ValueError: If the input 'image' is not an RGB image (3-channel).
    """
    # Check if image is an RGB image
    if len(image.shape) != 3 or image.shape[2] != 3:
        raise ValueError("The input 'image' must be an XxYxC RGB image - with three color channels (C).")

    # Check if flip argument is a list with two elements
    if not isinstance(flip, list) or len(flip) != 2:
        raise ValueError("The 'flip' argument should be a list with two elements [x_flip, y_flip].")

    # Check if crop argument is a list with four elements
    if not isinstance(crop, list) or len(crop) != 4:
        raise ValueError("The 'crop' argument should be a list with four elements [x_min, x_max, y_min, y_max].")

    hsv_image = rgb2hsv(image)

    if mask_tissue:
        s_image_gaussian = gaussian(hsv_image[..., 1], sigma=mask_gaussian_blur)
        thresh = threshold_otsu(s_image_gaussian)
        s_image_gaussian_binary = s_image_gaussian > thresh
        s_image_gaussian_binary = ndimage.binary_fill_holes(s_image_gaussian_binary).astype(int)
        image_out = image * s_image_gaussian_binary[..., np.newaxis]
        hsv_image_out = hsv_image * s_image_gaussian_binary[..., np.newaxis]

        if not keep_black_background:
            image_out = np.where(
                image_out == np.array([[0, 0, 0]]),
                np.array([[255, 255, 255]]),
                image_out,
            )
            hsv_image_out = np.where(
                hsv_image_out == np.array([[0, 0, 0]]),
                np.array([[1, 1, 1]]),
                hsv_image_out,
            )
    else:
        image_out = image
        hsv_image_out = hsv_image

    def _transform(im, flip, crop):
        return gaussian(rotate(equalize_adapthist(im)[:: flip[0], :: flip[1]].copy(), rotation), gaussian_blur)[
            crop[0] : crop[1], crop[2] : crop[3]
        ]

    prepared_images = [_transform(im, flip, crop) for im in hsv_image_out.transpose(2, 0, 1)]
    prepared_images += [_transform(im, flip, crop) for im in image_out.transpose(2, 0, 1)]
    return prepared_images


def prepare_pseudoimage_for_feature_matching(
    image: np.ndarray,
    invert: bool = True,
    gaussian_blur: float = 0,
    flip: list = [1, 1],
) -> np.ndarray:
    """
    Prepare an image for feature matching by applying optional transformations.

    Args:
        image (np.ndarray): Input image to be prepared.
        invert (bool, options): The image values will be inverted (white/black background).
        gaussian_blur (float, optional): Standard deviation for Gaussian blurring. Default is 0 (no blur).
        flip (list, optional): List indicating whether to flip the image along the x and y axes. Default is [1, 1].

    Returns:
        list: A list containing the prepared image after applying transformations.
            - The list contains a single element, which is the processed image.

    Notes:
        - This function applies optional Gaussian blurring and flipping to the input image.
        - The 'gaussian_blur' parameter controls the amount of blurring applied to the image.
        - The 'flip' parameter specifies flipping along the x and y axes using a list of factors [x_flip, y_flip].
    """
    _image = image.copy() if not invert else ((-image) - ((-image).min()))
    return [gaussian(_image, gaussian_blur)[:: flip[0], :: flip[1]]]


def run_registration(
    in_coords: np.ndarray,
    total_counts: np.ndarray,
    puck_id: np.ndarray,
    staining_image: np.ndarray,
    args,
) -> (np.ndarray, np.ndarray, PairwiseAlignmentMetadata):
    """
    Perform registration of spatial transcriptomics (STS) data with a staining image.

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
    # Create metadata object
    metadata = PairwiseAlignmentMetadata(args)

    # Create output object
    out_coords_output_fine = np.zeros_like(in_coords)

    # STAGE 1: coarse alignment (low-resolution feature matching)
    logging.info(f"Coarse registration, {len(in_coords)} coordinates")

    # Filter coordinates to only include UMI > args.threshold_counts_coarse
    sts_coords = in_coords[total_counts > args.threshold_counts_coarse]

    staining_image_rescaled = staining_image[:: args.rescale_factor_coarse, :: args.rescale_factor_coarse]

    # Create a pseudoimage from the STS coordinates
    sts_pseudoimage = create_pseudoimage(sts_coords, args.pseudoimage_size_coarse, staining_image_rescaled.shape)

    # Run feature matching between modalities
    dst = prepare_pseudoimage_for_feature_matching(sts_pseudoimage["pseudoimage"])

    _flips = [[1, 1], [1, -1], [-1, 1], [-1, -1]]
    _rotations = [0] if args.feature_matcher in ROTATION_INVARIANT else [0, 90]

    max_keypoints = 0
    _best_flip = _flips[0]
    _best_rotation = _rotations[0]
    _best_mkpts0 = None
    _best_mkpts1 = None

    for (_flip_x, _flip_y), _rotation in product(_flips, _rotations):
        logging.info(f"Aligning with image modality at flip {_flip_x, _flip_y} and rotation {_rotation}")
        # Preparing image and pseudoimage modalities for the feature matching model
        src = prepare_image_for_feature_matching(
            staining_image_rescaled,
            flip=[_flip_x, _flip_y],
            rotation=_rotation,
            mask_tissue=args.mask_tissue,
            keep_black_background=args.keep_black_background,
            mask_gaussian_blur=args.tissue_masking_gaussian_sigma,
        )

        # Find matching keypoints between image and STS pseudoimage modalities
        mkpts0, mkpts1 = feature_matching.find_matches(src, dst, args.feature_matcher)

        logging.info(f"{len(mkpts0)} matches")

        # Run RANSAC to remove outliers
        _, inliers = ransac(
            (mkpts0, mkpts1),
            SimilarityTransform,
            min_samples=args.ransac_coarse_min_samples,
            residual_threshold=args.ransac_coarse_residual_threshold,
            max_trials=args.ransac_coarse_max_trials,
        )
        inliers = inliers > 0

        if len(mkpts0) > max_keypoints:
            max_keypoints = len(mkpts0)
            _best_flip = [_flip_x, _flip_y]
            _best_rotation = _rotation
            _best_mkpts0 = mkpts0[inliers.flatten()]
            _best_mkpts1 = mkpts1[inliers.flatten()]

        logging.info(f"{inliers.sum()} inliers (RANSAC)")

    # Retrieve the results for the best flip combination
    # Filter keypoints with selected inliers
    in_mkpts0 = _best_mkpts0
    in_mkpts1 = _best_mkpts1

    # Estimate the transform
    tform_points = estimate_transform("similarity", in_mkpts0, in_mkpts1)

    # Apply the transform
    sts_coords_to_transform = sts_pseudoimage["coords_rescaled"] * sts_pseudoimage["rescaling_factor"]
    sts_coords_transformed = np.dot(
        tform_points.params,
        np.concatenate(
            [
                sts_coords_to_transform[:, ::-1],
                np.ones((len(sts_coords_to_transform), 1)),
            ],
            axis=1,
        ).T,
    ).T
    sts_coords_transformed = sts_coords_transformed * args.rescale_factor_coarse

    # Filter out the sts locations that lie outside of the staining image
    _i_sts_coords_coarse_within_image_bounds = np.where(
        (sts_coords_transformed[:, 0] > 0)
        & (sts_coords_transformed[:, 0] < staining_image.shape[1])
        & (sts_coords_transformed[:, 1] > 0)
        & (sts_coords_transformed[:, 1] < staining_image.shape[0])
    )
    sts_coords_transformed = sts_coords_transformed[_i_sts_coords_coarse_within_image_bounds][:, :-1]

    # Transform 'all' coordinates, to accomodate the fine transformation
    sts_coords_fine = in_coords.copy()
    sts_coords_fine -= sts_pseudoimage["offset_factor"]
    sts_coords_fine = (
        (sts_coords_fine / sts_pseudoimage["rescale_factor"]) * sts_pseudoimage["scale"]
    ) * sts_pseudoimage["rescaling_factor"]

    sts_coords_fine = np.dot(
        tform_points.params,
        np.concatenate(
            [
                sts_coords_fine[:, ::-1],
                np.ones((len(sts_coords_fine), 1)),
            ],
            axis=1,
        ).T,
    ).T
    sts_coords_fine = sts_coords_fine * args.rescale_factor_coarse
    sts_coords_fine = sts_coords_fine[:, :-1]
    out_coords_output_coarse = sts_coords_fine.copy()

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

    if args.only_coarse:
        return (
            out_coords_output_coarse,
            None,
            rotate(staining_image[:: _best_flip[0], :: _best_flip[1]], _best_rotation),
            metadata,
        )

    # STAGE 2: fine registration per tile
    logging.info(f"Fine registration with {_best_flip} flip and {_best_rotation} rotation")

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
        )

        # Axis limits to crop both modalities to tile region
        _t_sts_coords_to_transform = _t_sts_pseudoimage["coords_rescaled"] * _t_sts_pseudoimage["rescaling_factor"]

        min_lim, max_lim = _t_sts_coords_to_transform[_t_valid_coords].min(axis=0).astype(
            int
        ), _t_sts_coords_to_transform[_t_valid_coords].max(axis=0).astype(int)
        x_min, y_min = min_lim
        x_max, y_max = max_lim

        # Preparing image and pseudoimage modalities for feature detection (imaging modality has optimal flip)
        src = prepare_image_for_feature_matching(
            staining_image_rescaled,
            flip=_best_flip,
            rotation=_best_rotation,
            gaussian_blur=args.fine_registration_gaussian_sigma,
            crop=[x_min, x_max, y_min, y_max],
            mask_tissue=args.mask_tissue,
            keep_black_background=args.keep_black_background,
            mask_gaussian_blur=args.tissue_masking_gaussian_sigma,
        )
        dst = prepare_pseudoimage_for_feature_matching(
            _t_sts_pseudoimage["pseudoimage"][x_min:x_max, y_min:y_max],
            gaussian_blur=args.fine_registration_gaussian_sigma,
        )

        # Finding matches between modalities
        _t_mkpts0, _t_mkpts1 = feature_matching.find_matches(src, dst, args.feature_matcher)
        logging.info(f"{len(_t_mkpts0)} matches*")

        if len(_t_mkpts0) < args.fine_min_matches:
            logging.warn(f"Less than {args.fine_min_matches} were found. Skipping fine alignment at {tile_code}.")
            out_coords_output_fine[_t_puck_id] = sts_coords_fine[_t_puck_id]
            continue

        _, inliers = ransac(
            (_t_mkpts0, _t_mkpts1),
            SimilarityTransform,
            min_samples=args.ransac_fine_min_samples,
            residual_threshold=args.ransac_fine_residual_threshold,
            max_trials=args.ransac_fine_max_trials,
        )
        inliers = inliers > 0
        logging.info(f"{inliers.sum()} inliers (RANSAC)")

        # Compute similarity matrix and compute point transformation
        _t_tform_points = estimate_transform("similarity", _t_mkpts0[inliers.flatten()], _t_mkpts1[inliers.flatten()])

        # Apply the same transformation to the tiles
        _t_sts_coords_fine_to_transform = sts_coords_fine[_t_puck_id] / args.rescale_factor_fine

        # Check if transform within the acceptable bounds
        if (
            (_t_tform_points.rotation > np.pi / 4)
            or (_t_tform_points.scale > 2 or _t_tform_points.scale < 0.5)
            or (_t_tform_points.translation.max() > _t_sts_coords_fine_to_transform.max(axis=0).max())
        ):
            logging.warn(f"Transformation matrix out of bounds for tile {tile_code}")
            out_coords_output_fine[_t_puck_id] = sts_coords_fine[_t_puck_id]
            continue

        # If the previous filter passes, apply the transformation
        _t_sts_coords_fine_to_transform = (_t_sts_coords_fine_to_transform - np.array([[y_min, x_min]]))[:, ::-1]
        _t_sts_coords_fine_transformed = np.dot(
            _t_tform_points.params,
            np.concatenate(
                [
                    _t_sts_coords_fine_to_transform[:, ::-1],
                    np.ones((len(_t_sts_coords_fine_to_transform), 1)),
                ],
                axis=1,
            ).T,
        ).T

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
        rotate(staining_image[:: _best_flip[0], :: _best_flip[1]], _best_rotation),
        metadata,
    )


def _run_pairwise_aligner(args):
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
    sts_aligned_coarse, sts_aligned_fine, staining_image_aligned, metadata = run_registration(
        sts["obsm/spatial"],
        sts["obs/total_counts"],
        sts["obs/puck_id"],
        staining_image,
        args,
    )

    # Saving the metadata (for QC)
    if args.metadata_out != "":
        save_pickle(metadata, args.metadata_out)

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
    args = get_pairwise_aligner_parser().parse_args()
    with threadpool_limits(limits=args.n_threads):
        _run_pairwise_aligner(args)
