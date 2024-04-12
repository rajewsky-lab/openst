import json
import logging
import os

import anndata as ad
import numpy as np
from anndata import AnnData, read_h5ad
from PIL import Image
# from scipy.ndimage import binary_fill_holes
# from skimage.color import rgb2hsv
# from skimage.filters import gaussian, threshold_otsu
from skimage.transform import estimate_transform, resize, warp
from tqdm import tqdm

from openst.utils.file import check_directory_exists, check_file_exists


def transform_coordinates_h5_from_stim(adata: AnnData, transform_model: np.ndarray, rescale: float = 1):
    """
    Transform spatial coordinates in an AnnData object using an affine transformation matrix.

    Args:
        adata (AnnData): Annotated data object containing spatial coordinates to be transformed.
        transform_model (np.ndarray): Affine transformation matrix to apply to the spatial coordinates.
        rescale (float, optional): Scaling factor to rescale the transformed coordinates. Defaults to 1.

    Returns:
        AnnData: Annotated data object with transformed spatial coordinates.

    Notes:
        - This function transforms spatial coordinates using an affine transformation matrix.
        - The 'rescale' parameter allows rescaling of the transformed coordinates.
    """
    adata.obsm["spatial_transform"] = (
        np.dot(
            transform_model,
            np.concatenate([adata.obsm["spatial"], np.ones((len(adata.obsm["spatial"]), 1))], axis=1).T,
        ).T
        / rescale
    )
    adata.uns["STIM_aligned"] = {"affine_matrix": transform_model, "rescale_factor": rescale}

    return adata


def estimate_transform_for_image_from_stim(
    image: np.ndarray, section: AnnData, transform_model: np.ndarray, rescale: float = 1
):
    """
    Estimate transformed coordinates for an image based on an AnnData object and an affine transformation.

    Args:
        image (np.ndarray): Input image for which transformed coordinates need to be estimated.
        section (AnnData): AnnData object representing a section containing spatial information.
        transform_model (np.ndarray): Affine transformation matrix used for estimation.
        rescale (float, optional): Scaling factor used for rescaling coordinates. Defaults to 1.

    Returns:
        np.ndarray: Transformed coordinates for the image.
        np.ndarray: Meshgrid coordinates used for estimation.

    Notes:
        - This function estimates transformed coordinates for an image based
          on an AnnData object and a transformation matrix.
        - The 'rescale' parameter allows rescaling of the estimated coordinates.
    """
    x = np.arange(image.shape[1])
    y = np.arange(image.shape[0])
    xv, yv = np.meshgrid(x, y)

    mgrid = np.concatenate([xv[..., np.newaxis], yv[..., np.newaxis]], axis=2).reshape(-1, 2)[..., ::-1]
    mgrid = mgrid[[0, (mgrid.size * 0.25).astype(int), (mgrid.size * 0.75).astype(int), -1]]

    coords_orig_im = (mgrid - section["spatial"].min(axis=0)).astype(int) * rescale
    coords_transform_image = (
        np.dot(transform_model, np.concatenate([coords_orig_im, np.ones((len(coords_orig_im), 1))], axis=1).T).T
        / rescale
    )

    return coords_transform_image, mgrid


def transform_images_from_stim_results(
    models_sift: list, sample_aligned: list, images_sections: list, downsample_image: int = 1
):
    """
    Transform images based on provided SIFT models, AnnData objects, and image sections.

    Args:
        models_sift (list): List of SIFT models for image transformation.
        sample_aligned (list): List of AnnData objects representing aligned samples.
        images_sections (list): List of images corresponding to the sections.
        downsample_image (int, optional): Downsampling factor for image processing. Defaults to 1.

    Returns:
        list: List of transformed images.

    Notes:
        - This function transforms images using SIFT models, AnnData objects, and image sections.
        - The 'downsample_image' parameter controls downsampling of the images for processing.
    """
    if not isinstance(downsample_image, int) or downsample_image < 1:
        raise ValueError("'downsample_image' needs to be an integer larger than 1")

    if len(models_sift) != len(sample_aligned) and len(sample_aligned) != len(images_sections):
        raise ValueError("The objects 'models_sift', 'sample_aligned' and 'images_sections' need to have same length")

    # Create image output objects
    section_affine_transformed_mgrid = []
    section_affine_transformed_coords_image = []
    section_affine_transformed_image = []

    # Calculate affine transformation for images
    for model_sift, section, image in tqdm(models_sift, sample_aligned, images_sections):
        coords_transform_image, mgrid = estimate_transform_for_image_from_stim(image, model_sift, section)
        section_affine_transformed_mgrid.append(mgrid)
        section_affine_transformed_coords_image.append(coords_transform_image)

    min_coord, max_coord = np.concatenate(section_affine_transformed_coords_image, axis=0).min(axis=0), np.concatenate(
        section_affine_transformed_coords_image, axis=0
    ).max(axis=0)

    # Apply the affine transformation to the images
    for section, image, coord_image, mgrid in tqdm(
        sample_aligned, images_sections, section_affine_transformed_coords_image, section_affine_transformed_mgrid
    ):
        image_to_transform = resize(
            image,
            image[:: args.downsample_image, :: args.downsample_image].shape,
            anti_aliasing=True,
            preserve_range=True,
        ).astype(np.uint8)
        tform = estimate_transform(
            "affine",
            (coord_image[..., :2] - min_coord[:2][np.newaxis]) / args.downsample_image,
            mgrid / args.downsample_image,
        )
        images_sections_transformed_tformed = warp(
            image_to_transform.transpose(1, 0, 2),
            tform,
            mode="constant",
            preserve_range=True,
            output_shape=((max_coord / args.downsample_image).astype(int)),
            order=1,
        )
        section_affine_transformed_image.append(images_sections_transformed_tformed)

    return section_affine_transformed_image


def _run_from_3d_registration(args):
    """
    _run_from_3d_registration
    """

    # Check input and output data
    check_file_exists(args.in_adata_aligned)

    if not check_directory_exists(args.h5_out):
        raise FileNotFoundError("Parent directory for --output-mask does not exist")

    Image.MAX_IMAGE_PIXELS = args.max_image_pixels

    # Loading input
    # Loading SIFT model from STIM
    samples_model_sift = []
    for sample in tqdm(args.n5_dirs):
        with open(os.path.join(sample, "attributes.json"), "r") as f:
            n5_json_object = json.load(f)
            model_sift = n5_json_object["model_sift"]
            model_sift = np.array(
                [
                    [model_sift[0], model_sift[1], model_sift[2]],
                    [model_sift[3], model_sift[4], model_sift[5]],
                    [0, 0, 1],
                ]
            )
            samples_model_sift.append(model_sift)

    # Loading spatial transcriptomics data
    samples_adata = []
    for sample in tqdm(args.h5_files):
        samples_adata.append(read_h5ad(sample))

    # Loading image data, if specified
    if args.images == "":
        images_sections = None
        if not check_directory_exists(args.image_out):
            raise FileNotFoundError("Parent directory for --output-mask does not exist")

    if args.images_in_adata:
        if isinstance(args.images, list) and len(args.images) > 1:
            raise ValueError("--images must be of length 1 when --images-in-adata is specified.")

        images_sections = [sample[args.images] for sample in samples_adata]
    else:
        for image in args.images:
            check_file_exists(image)

        images_sections = [np.array(Image.open(image)) for image in args.images]

    # Perform coordinate transformation of spatial transcriptomics data
    sample_aligned = []
    for model_sift, sample, name in tqdm(zip(samples_model_sift, samples_adata, args.h5_files)):
        sample.obsm["spatial"] -= sample.obsm["spatial"].min(axis=0)

        if args.rescale != 1:
            sample.obsm["spatial"] = sample.obsm["spatial"] * args.rescale

        sample = transform_coordinates_h5_from_stim(sample, model_sift, args.rescale)
        sample.obs["n_section"] = name
        sample.obs_names = sample.obs["n_section"].astype(str) + ":" + sample.obs_names.astype(str)
        sample_aligned.append(sample)

    adata_aligned = ad.concat(sample_aligned, merge=args.merge_output, join=args.join_output)
    adata_aligned.obs_names_make_unique()

    if images_sections is not None:
        # Register images to STIM coordinates
        section_affine_transformed_image = transform_images_from_stim_results(
            samples_model_sift, sample_aligned, images_sections, downsample_image=args.downsample_image
        )

        # Save image data according to specified arguments
        if not args.save_images_separately:
            adata_aligned.uns["STIM_aligned"]["staining_images"] = {}
        for image, name in zip(section_affine_transformed_image, args.h5_files):
            if args.save_images_separately:
                Image.fromarray(image.astype(np.uint8)).save(os.path.join(args.image_out, f"{name}_aligned.tif"))
            else:
                adata_aligned.uns["STIM_aligned"]["staining_images"][name] = image

    # Save final adata object
    adata_aligned.write(args.h5_out)


if __name__ == "__main__":
    from openst.cli import get_from_3d_registration_parser
    args = get_from_3d_registration_parser().parse_args()
    _run_from_3d_registration(args)
