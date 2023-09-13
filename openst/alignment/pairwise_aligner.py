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
from itertools import product
import logging

import numpy as np
from anndata import read_h5ad
from PIL import Image
from scipy import ndimage
from skimage.color import rgb2hsv, rgb2gray
from skimage.exposure import equalize_adapthist
from skimage.filters import gaussian, threshold_otsu
from skimage.transform import estimate_transform, rotate
from threadpoolctl import threadpool_limits

from openst.alignment import feature_matching, fiducial_detection
from openst.metadata.classes.pairwise_alignment import (AlignmentResult,
                                       PairwiseAlignmentMetadata)
from openst.alignment.transformation import apply_transform
from openst.utils.file import (check_adata_structure, check_directory_exists,
                               check_file_exists, load_properties_from_adata,)
from openst.utils.pseudoimage import create_pseudoimage


def get_pairwise_aligner_parser():
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
        default=50,
        help="Times RANSAC will run (x1000 iterations) during coarse registration",
    )
    parser.add_argument(
        "--ransac-fine-min-samples",
        type=int,
        default=3,
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
        default=50,
        help="Times RANSAC will run (x1000 iterations) during fine registration",
    )
    parser.add_argument(
        "--max-image-pixels",
        type=int,
        default=933120000,
        help="Upper bound for number of pixels in the images (prevents exception when opening very large images)",
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
    parser.add_argument(
        "--fiducial-model",
        type=str,
        default="",
        help="Path to a object detection model (YOLO) to detect fiducial markers",
    )
    parser.add_argument(
        "--genes-coarse",
        nargs="+",
        type=str,
        default=None,
        help="Genes used for plotting the pseudoimage during the coarse alignment phase.",
    )
    parser.add_argument(
        "--genes-fine",
        nargs="+",
        type=str,
        default=None,
        help="Genes used for plotting the pseudoimage during the fine alignment phase.",
    )
    parser.add_argument(
        "--device",
        type=str,
        default='cpu',
        choices=['cpu', 'cuda'],
        help="Device used to run feature matching model. Can be ['cpu', 'cuda']",
    )
    return parser


def setup_pairwise_aligner_parser(parent_parser):
    """setup_pairwise_aligner_parser"""
    parser = parent_parser.add_parser(
        "pairwise_aligner",
        help="openst pairwise alignment of two-dimensional spatial transcriptomics and imaging data",
        parents=[get_pairwise_aligner_parser()],
    )
    parser.set_defaults(func=_run_pairwise_aligner)

    return parser


def transform_image(im, flip: list = None, crop: list = None, rotation: int = None):
    _im = im
    if flip is not None:
        _im = _im[:: flip[0], :: flip[1]]
    if rotation is not None:
        _im = rotate(_im, rotation)
    if crop is not None:
        _im = _im[crop[0] : crop[1], crop[2] : crop[3]]
    
    return _im


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

    image = transform_image(image, flip, crop, rotation)
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

    prepared_images = [gaussian((equalize_adapthist(rgb2gray(image_out))), gaussian_blur)]
    prepared_images += [gaussian((equalize_adapthist(im)), gaussian_blur) for im in hsv_image_out.transpose(2, 0, 1)]
    prepared_images += [gaussian((equalize_adapthist(im)), gaussian_blur) for im in image_out.transpose(2, 0, 1)]
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
    return [gaussian(equalize_adapthist(_image), gaussian_blur)[:: flip[0], :: flip[1]]]


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
    # Create output objects
    out_coords_output_fine = np.zeros_like(in_coords)
    metadata = PairwiseAlignmentMetadata(args)

    logging.info(f"Coarse registration, {len(in_coords)} coordinates")

    # Preparing images and preprocessing routines
    sts_coords = in_coords[total_counts > args.threshold_counts_coarse]
    staining_image_rescaled = staining_image[:: args.rescale_factor_coarse, :: args.rescale_factor_coarse]
    sts_pseudoimage = create_pseudoimage(sts_coords, args.pseudoimage_size_coarse, staining_image_rescaled.shape)

    def src_augmenter(x, flip, rotation):
        return prepare_image_for_feature_matching(
            image=x,
            flip=flip,
            rotation=rotation,
            mask_tissue=args.mask_tissue,
            keep_black_background=args.keep_black_background,
            mask_gaussian_blur=args.tissue_masking_gaussian_sigma,
        )

    dst = prepare_pseudoimage_for_feature_matching(sts_pseudoimage["pseudoimage"])
    src = staining_image_rescaled

    # Feature matching
    in_mkpts0, in_mkpts1, _best_flip, _best_rotation = feature_matching.match_images(
        src,
        dst,
        feature_matcher=args.feature_matcher,
        src_augmenter=src_augmenter,
        ransac_min_samples=args.ransac_coarse_min_samples,
        ransac_residual_threshold=args.ransac_coarse_residual_threshold,
        ransac_max_trials=args.ransac_coarse_max_trials,
        device=args.device
    )

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
            rotate(staining_image[:: _best_flip[0], :: _best_flip[1]], _best_rotation),
            metadata,
        )

    # STAGE 2: fine registration per tile
    logging.info(f"Fine registration with {_best_flip} flip and {_best_rotation} rotation")

    # Collect tile identifiers
    tile_codes = np.unique(puck_id.codes)

    # Apply scaling to input image again, for fine registration
    staining_image_rescaled = staining_image[:: args.rescale_factor_fine, :: args.rescale_factor_fine]
    src = staining_image_rescaled
    src = transform_image(staining_image_rescaled, _best_flip, None, _best_rotation)

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
            values=None
        )

        _t_counts = total_counts[(total_counts > args.threshold_counts_coarse)][_i_sts_coords_coarse_within_image_bounds][_t_valid_coords]

        _t_sts_pseudoimage_counts = create_pseudoimage(
            sts_coords_transformed[:, ::-1],  # we need to flip these coordinates
            args.pseudoimage_size_fine,
            staining_image_rescaled.shape,
            _t_valid_coords,
            recenter=False,
            rescale=True,
            values=_t_counts.astype(int)
        )

        # Axis limits to crop both modalities to tile region
        _t_sts_coords_to_transform = _t_sts_pseudoimage["coords_rescaled"] * _t_sts_pseudoimage["rescaling_factor"]

        min_lim, max_lim = _t_sts_coords_to_transform[_t_valid_coords].min(axis=0).astype(
            int
        ), _t_sts_coords_to_transform[_t_valid_coords].max(axis=0).astype(int)
        x_min, y_min = min_lim
        x_max, y_max = max_lim

        # Preparing image and pseudoimage modalities for feature detection (imaging modality has optimal flip)
        def src_preprocessor(x, flip, rotation):
            return prepare_image_for_feature_matching(
                image=x,
                gaussian_blur=args.fine_registration_gaussian_sigma,
                crop=[x_min, x_max, y_min, y_max],
                mask_tissue=args.mask_tissue,
                keep_black_background=args.keep_black_background,
                mask_gaussian_blur=args.tissue_masking_gaussian_sigma,
            )

        dst = []
        for pseudoimage, invert in product([_t_sts_pseudoimage, _t_sts_pseudoimage_counts], [True, False]):
            dst = prepare_pseudoimage_for_feature_matching(
                pseudoimage["pseudoimage"][x_min:x_max, y_min:y_max],
                gaussian_blur=args.fine_registration_gaussian_sigma,
                invert=invert
            )

        # Finding matches between modalities
        _t_mkpts0, _t_mkpts1, _, _ = feature_matching.match_images(
            src,
            dst,
            flips=[_best_flip],
            rotations=[_best_rotation],
            src_augmenter=src_preprocessor,
            ransac_min_samples=args.ransac_fine_min_samples,
            ransac_residual_threshold=args.ransac_fine_residual_threshold,
            ransac_max_trials=args.ransac_fine_max_trials,
            device=args.device
        )

        # Detect fiducial markers with the YOLO model
        if args.fiducial_model != "" and check_file_exists(args.fiducial_model, exception=False):
            logging.info("Running fiducial marker refinement")
            fiducial_points_src = fiducial_detection.find_fiducial(
                transform_image(src, _best_flip, [x_min, x_max, y_min, y_max], _best_rotation), args.fiducial_model
            )
            fiducial_points_dst = fiducial_detection.find_fiducial(dst, args.fiducial_model)

            _t_mkpts0_fiducial, _t_mkpts1_fiducial = fiducial_detection.correspondences_fiducials(
                fiducial_points_src, fiducial_points_dst, distance_threshold=50
            )

            logging.info(f"{len(_t_mkpts0_fiducial)} fiducial matches")
            _t_mkpts0 = np.concatenate([_t_mkpts0_fiducial, _t_mkpts0], axis=1)
            _t_mkpts1 = np.concatenate([_t_mkpts1_fiducial, _t_mkpts1], axis=1)

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
        src.astype(np.uint8),
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
    args = get_pairwise_aligner_parser().parse_args()
    with threadpool_limits(limits=args.n_threads):
        _run_pairwise_aligner(args)
