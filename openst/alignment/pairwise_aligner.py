"""
Automatic Pairwise Alignment of Spatial Transcriptomics and Imaging Data (Open-ST)

Author: Daniel León-Periñán @ N.Rajewsky Lab (BIMSB)
Date: April 15, 2024
Version: 0.2.0

Description:
This script performs automatic pairwise alignment between spatial transcriptomics data and imaging data.
It aligns the spatial transcriptomics data onto the imaging data to enable spatial comparison and analysis.
"""

import logging
from itertools import product

import cv2
import numpy as np
import h5py
from skimage.color import rgb2gray, rgb2hsv
from skimage.exposure import equalize_adapthist
from skimage.filters import gaussian
from skimage.transform import estimate_transform, rescale, rotate
from threadpoolctl import threadpool_limits

from openst.alignment import feature_matching, fiducial_detection
from openst.alignment.transformation import apply_transform
from openst.metadata.classes.pairwise_alignment import (
    AlignmentResult, PairwiseAlignmentMetadata)
from openst.utils.file import (check_adata_structure, check_directory_exists,
                               check_file_exists, load_properties_from_adata)
from openst.utils.pimage import mask_tissue as p_mask_tissue
from openst.utils.pimage import is_grayscale
from openst.utils.pseudoimage import create_paired_pseudoimage


def transform_image(im, flip: list = None, crop: list = None, rotation: int = None):
    _dtype = im.dtype
    _im = im
    if flip is not None:
        _im = _im[:: flip[0], :: flip[1]]
    if rotation is not None:
        _im = rotate(_im, rotation, clip=True, preserve_range=True, resize=True)
    if crop is not None:
        _im = _im[crop[0] : crop[1], crop[2] : crop[3]]

    return _im.astype(_dtype)


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
    
    if mask_tissue:
        image_out, hsv_image_out = p_mask_tissue(image,
                                               keep_black_background,
                                               mask_gaussian_blur,
                                               return_hsv=True)
    else:
        image_out = image
        hsv_image_out = rgb2hsv(image)

    prepared_images = [gaussian((equalize_adapthist(rgb2gray(image_out))), gaussian_blur)]
    prepared_images += [gaussian((equalize_adapthist(im)), gaussian_blur) for im in hsv_image_out.transpose(2, 0, 1)]
    prepared_images += [gaussian((equalize_adapthist(im)), gaussian_blur) for im in image_out.transpose(2, 0, 1)]
    return prepared_images


def prepare_image_for_feature_matching_grayscale(
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
    # Check if flip argument is a list with two elements
    if not isinstance(flip, list) or len(flip) != 2:
        raise ValueError("The 'flip' argument should be a list with two elements [x_flip, y_flip].")

    # Check if crop argument is a list with four elements
    if not isinstance(crop, list) or len(crop) != 4:
        raise ValueError("The 'crop' argument should be a list with four elements [x_min, x_max, y_min, y_max].")

    image = transform_image(image, flip, crop, rotation)
    image_out = image

    prepared_images = [gaussian((equalize_adapthist(image_out)), gaussian_blur)]
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
    tile_id: np.ndarray,
    staining_image: np.ndarray,
    args,
) -> (np.ndarray, np.ndarray, np.ndarray, PairwiseAlignmentMetadata): 
    """
    Perform registration of spatial transcriptomics (STS) data with a staining image.

    Args:
        in_coords (np.ndarray): Input STS coordinates.
        total_counts (np.ndarray): Total UMI counts for each STS coordinate.
        tile_id: Identifier for each STS coordinate. During the fine registration,
                 this 'tile_id' is used to aggregate the coordinates into buckets that
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
    
    logging.info(f"Rescaling input image for coarse registration")
    # rescaling image
    factors = np.array([args.rescale_factor_coarse, args.rescale_factor_coarse])
    _ker = np.maximum(0, (factors - 1) / 2).astype(int)
    src = cv2.resize(cv2.blur(staining_image, _ker),
                     (np.array(staining_image.shape)[[0, 1]]/args.rescale_factor_coarse).astype(int)[::-1],
                     interpolation=cv2.INTER_NEAREST)
    # staining_image_rescaled = rescale(
    #     staining_image, 1 / args.rescale_factor_coarse, preserve_range=True, anti_aliasing=True, channel_axis=-1
    # ).astype(np.uint8)
    # src = staining_image_rescaled

    _fn_prepare_image_for_feature_matching = prepare_image_for_feature_matching

    if is_grayscale(src):
        _fn_prepare_image_for_feature_matching = prepare_image_for_feature_matching_grayscale

    def src_augmenter(x, flip, rotation):
        return _fn_prepare_image_for_feature_matching(
            image=x,
            flip=flip,
            rotation=rotation,
            mask_tissue=args.mask_tissue,
            keep_black_background=args.keep_black_background,
            mask_gaussian_blur=args.tissue_masking_gaussian_sigma,
        )

    sts_coords = in_coords[total_counts > args.threshold_counts_coarse]
    sts_pseudoimage = create_paired_pseudoimage(sts_coords, args.pseudoimage_size_coarse, src.shape, resize_method='cv2')
    dst = prepare_pseudoimage_for_feature_matching(sts_pseudoimage["pseudoimage"])

    # Feature matching
    in_mkpts0, in_mkpts1, _best_flip, _best_rotation = feature_matching.match_images(
        src,
        dst,
        feature_matcher=args.feature_matcher,
        src_augmenter=src_augmenter,
        ransac_min_samples=args.ransac_coarse_min_samples,
        ransac_residual_threshold=args.ransac_coarse_residual_threshold,
        ransac_max_trials=args.ransac_coarse_max_trials,
        device=args.device,
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
    sts_coords_transformed = sts_coords_transformed[_i_sts_coords_coarse_within_image_bounds][:, [0, 1]]

    # Apply transform to all coordinates & retransform back
    sts_coords_coarse = in_coords.copy()
    sts_coords_coarse -= sts_pseudoimage["offset_factor"]
    sts_coords_coarse = (
        (sts_coords_coarse / sts_pseudoimage["rescale_factor"]) * sts_pseudoimage["scale"]
    ) * sts_pseudoimage["rescaling_factor"]

    sts_coords_coarse = apply_transform(sts_coords_coarse, tform_points, check_bounds=False)
    sts_coords_coarse = sts_coords_coarse * args.rescale_factor_coarse
    sts_coords_coarse = sts_coords_coarse[:, [0, 1]]
    out_coords_output_coarse = sts_coords_coarse.copy()

    # Saving alignment results here
    # TODO: check order of keypoints (in all functions throughout package)
    _align_result = AlignmentResult(
        name="coarse_alignment_whole_section",
        im_0=transform_image(src, _best_flip, [0, None, 0, None], _best_rotation),
        im_1=sts_pseudoimage["pseudoimage"],
        transformation_matrix=tform_points.params,
        ransac_results=None,
        sift_results=None,
        keypoints0=in_mkpts1,
        keypoints1=in_mkpts0,
    )
    metadata.add_alignment_result(_align_result)

    # Finish here if only coarse registration was selected to run
    if args.only_coarse:
        return (
            out_coords_output_coarse,
            None,
            transform_image(staining_image, _best_flip, None, _best_rotation).astype(np.uint8),
            metadata,
        )

    # STAGE 2: fine registration per tile
    logging.info(f"Fine registration with {_best_flip} flip and {_best_rotation} rotation")

    # Collect tile identifiers
    tile_codes = np.unique(tile_id.codes)

    # Apply scaling to input image again, for fine registration
    staining_image_rescaled = rescale(
        staining_image, 1 / args.rescale_factor_fine, preserve_range=True, anti_aliasing=True, channel_axis=-1
    ).astype(np.uint8)
    src = staining_image_rescaled.astype(np.uint8)
    src = transform_image(staining_image_rescaled, _best_flip, None, _best_rotation)

    for tile_code in tile_codes:
        # Create a pseudoimage
        _t_tile_id = np.isin(tile_id.codes, tile_code)
        _t_valid_coords = np.isin(
            tile_id[(total_counts > args.threshold_counts_coarse)].codes[_i_sts_coords_coarse_within_image_bounds],
            tile_code,
        )

        logging.info(f"Registering tile {tile_code} with {_t_valid_coords.sum()} coordinates")

        _t_sts_pseudoimage = create_paired_pseudoimage(
            sts_coords_transformed[:, ::-1],  # we flip these coordinates
            args.pseudoimage_size_fine,
            staining_image_rescaled.shape,
            _t_valid_coords,
            recenter=False,
            rescale=True,
            values=None,
        )

        _t_counts = total_counts[(total_counts > args.threshold_counts_coarse)][
            _i_sts_coords_coarse_within_image_bounds
        ][_t_valid_coords]

        _t_sts_pseudoimage_counts = create_paired_pseudoimage(
            sts_coords_transformed[:, ::-1],  # we flip these coordinates
            args.pseudoimage_size_fine,
            staining_image_rescaled.shape,
            _t_valid_coords,
            recenter=False,
            rescale=True,
            values=_t_counts,
        )

        # Axis limits to crop both modalities to tile region
        _t_sts_coords_to_transform = _t_sts_pseudoimage["coords_rescaled"] * _t_sts_pseudoimage["rescaling_factor"]

        min_lim, max_lim = _t_sts_coords_to_transform[_t_valid_coords].min(axis=0).astype(
            int
        ), _t_sts_coords_to_transform[_t_valid_coords].max(axis=0).astype(int)
        x_min, y_min = min_lim
        x_max, y_max = max_lim

        _fn_prepare_image_for_feature_matching = prepare_image_for_feature_matching

        if is_grayscale(src):
            _fn_prepare_image_for_feature_matching = prepare_image_for_feature_matching_grayscale

        # Preparing image and pseudoimage modalities for feature detection (imaging modality has optimal flip)
        def src_preprocessor(x, flip, rotation):
            return _fn_prepare_image_for_feature_matching(
                image=x,
                gaussian_blur=args.fine_registration_gaussian_sigma,
                crop=[x_min, x_max, y_min, y_max],
                mask_tissue=args.mask_tissue,
                keep_black_background=args.keep_black_background,
                mask_gaussian_blur=args.tissue_masking_gaussian_sigma,
            )

        dst = []
        for pseudoimage, invert in product(
            [_t_sts_pseudoimage["pseudoimage"], _t_sts_pseudoimage_counts["pseudoimage"].astype(int)], [False, True]
        ):
            dst += prepare_pseudoimage_for_feature_matching(
                pseudoimage[x_min:x_max, y_min:y_max],
                gaussian_blur=args.fine_registration_gaussian_sigma,
                invert=invert,
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
            device=args.device,
        )

        # Detect fiducial markers with the YOLO model
        if args.fiducial_model != "" and check_file_exists(args.fiducial_model, exception=False):
            raise NotImplementedError("Fiducial marker refinement is not implemented yet")
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

        # Apply the same transformation to the tiles
        _t_sts_coords_fine_to_transform = sts_coords_coarse[_t_tile_id] / args.rescale_factor_fine
        _t_sts_coords_fine_to_transform = (_t_sts_coords_fine_to_transform - np.array([[y_min, x_min]]))[:, ::-1]

        # Compute similarity matrix and compute point transformation
        if len(_t_mkpts0) > args.fine_min_matches:
            _t_tform_points = estimate_transform("similarity", _t_mkpts0, _t_mkpts1)

            _t_sts_coords_fine_transformed = apply_transform(
                _t_sts_coords_fine_to_transform, _t_tform_points, check_bounds=True
            )[:, :-1]

            _tform_params = _t_tform_points.params
        else:
            logging.warning(f"There were not enough matching points ({len(_t_mkpts0)} out of selected {args.fine_min_matches})")
            _t_sts_coords_fine_transformed = _t_sts_coords_fine_to_transform[:, ::-1]
            _tform_params = None
            

        # Rescale points to original HE dimensions
        _t_sts_coords_fine_transformed = _t_sts_coords_fine_transformed + np.array([[y_min, x_min]])
        _t_sts_coords_fine_transformed = _t_sts_coords_fine_transformed * args.rescale_factor_fine

        out_coords_output_fine[_t_tile_id] = _t_sts_coords_fine_transformed[:, ::-1]

        # Saving alignment results here (only when passed)
        # TODO: check order of keypoints (in all functions throughout package)
        _align_result = AlignmentResult(
            name=f"fine_alignment_tile_{tile_code}",
            im_0=src[x_min:x_max, y_min:y_max],
            im_1=_t_sts_pseudoimage["pseudoimage"][x_min:x_max, y_min:y_max],
            transformation_matrix=_tform_params,
            ransac_results=None,
            sift_results=None,
            keypoints0=_t_mkpts1,
            keypoints1=_t_mkpts0,
        )
        metadata.add_alignment_result(_align_result)

    return (
        out_coords_output_coarse,
        out_coords_output_fine,
        transform_image(staining_image, _best_flip, None, _best_rotation).astype(np.uint8),
        metadata,
    )


def run_pairwise_aligner(args):
    # Check input and output data
    check_file_exists(args.h5_in)
    check_adata_structure(args.h5_in)

    if args.fiducial_model != "" and check_file_exists(args.fiducial_model, exception=False):
        raise NotImplementedError("Fiducial marker refinement is not implemented yet")

    if args.metadata != "" and not check_directory_exists(args.metadata):
        raise FileNotFoundError("Parent directory for the metadata does not exist")

    # Loading the spatial transcriptomics data
    sts = load_properties_from_adata(args.h5_in, properties=["obsm/spatial", "obs/total_counts", "obs/tile_id"])

    # Loading image data
    with h5py.File(args.h5_in, 'r') as adata:
        staining_image = adata[args.image_in][:]

    # Running registration
    sts_aligned_coarse, sts_aligned_fine, staining_image_aligned, metadata = run_registration(
        sts["obsm/spatial"],
        sts["obs/total_counts"],
        sts["obs/tile_id"],
        staining_image,
        args,
    )

    # Saving the metadata (for QC)
    if args.metadata != "":
        metadata.render()
        metadata.save_json(args.metadata)

    logging.info(f"Updating {args.h5_in} in place")
    with h5py.File(args.h5_in, 'r+') as adata:
        # TODO: to this with a function, instead
        # TODO: do not save image because we assume it is already there!
        if "uns/spatial_pairwise_aligned/staining_image" in adata:
            del adata["uns/spatial_pairwise_aligned/staining_image"]
        adata["uns/spatial_pairwise_aligned/staining_image"] = staining_image

        if "uns/spatial_pairwise_aligned/staining_image_transformed" in adata:
            del adata["uns/spatial_pairwise_aligned/staining_image_transformed"]
        adata["uns/spatial_pairwise_aligned/staining_image_transformed"] = staining_image_aligned

        if "obsm/spatial_pairwise_aligned_coarse" in adata:
            adata["obsm/spatial_pairwise_aligned_coarse"][:] = sts_aligned_coarse[..., ::-1]
        else:
            adata["obsm/spatial_pairwise_aligned_coarse"] = sts_aligned_coarse[..., ::-1]
        
        if sts_aligned_fine is not None:
            if "obsm/spatial_pairwise_aligned_fine" in adata:
                adata["obsm/spatial_pairwise_aligned_fine"][:] = sts_aligned_fine
            else:
                adata["obsm/spatial_pairwise_aligned_fine"] = sts_aligned_fine


def _run_pairwise_aligner(args):
    with threadpool_limits(limits=args.num_workers):
        run_pairwise_aligner(args)


if __name__ == "__main__":
    from openst.cli import get_pairwise_aligner_parser
    args = get_pairwise_aligner_parser().parse_args()
    _run_pairwise_aligner(args)
