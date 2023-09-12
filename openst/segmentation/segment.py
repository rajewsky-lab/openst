import argparse
import logging

import numpy as np
from cellpose import models
from PIL import Image
from skimage import measure
from skimage.segmentation import expand_labels
from tqdm import tqdm

from openst.utils.file import check_directory_exists, check_file_exists


def get_segment_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="segmentation of open-ST imaging data with cellpose",
        allow_abbrev=False,
        add_help=False,
    )
    parser.add_argument(
        "--image-in",
        type=str,
        required=True,
        help="Path to the input image.",
    )
    parser.add_argument(
        "--output-mask",
        type=str,
        required=True,
        help="Path to the output file where the mask will be stored",
    )
    parser.add_argument(
        "--model",
        type=str,
        default="",
        help=f"""cellpose model - either a path or a valid model name.
        Valid model names are {models.MODEL_NAMES}""",
    )
    parser.add_argument(
        "--flow-threshold",
        type=float,
        default=0.5,
        help="cellpose's 'flow_threshold' parameter",
    )
    parser.add_argument(
        "--cellprob-threshold",
        type=float,
        default=0,
        help="cellpose's 'cellprob_threshold' parameter",
    )
    parser.add_argument(
        "--diameter",
        type=float,
        default=20,
        help="cellpose's 'diameter' parameter",
    )
    parser.add_argument(
        "--tile-size",
        type=int,
        default=512,
        help="When prediction of the mask runs in separate tiles, this is the tile square size (in pixels)",
    )
    parser.add_argument(
        "--by-tiles",
        action="store_true",
        help="When specified, segmentation is computed at non-overlapping tiles of size '--tile-size'",
    )
    parser.add_argument(
        "--max-image-pixels",
        type=int,
        default=933120000,
        help="Upper bound for number of pixels in the images (prevents exception when opening very large images)",
    )
    parser.add_argument(
        "--gpu",
        action="store_true",
        help="When specified, a GPU will be used for prediction",
    )
    parser.add_argument(
        "--dilate-px",
        type=int,
        help="Pixels the outlines of the segmentation mask will be extended",
        required=False,
        default=0,
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


def setup_segment_parser(parent_parser):
    """setup_segment_parser"""
    parser = parent_parser.add_parser(
        "segment",
        help="segmentation of open-ST imaging data with cellpose",
        parents=[get_segment_parser()],
    )
    parser.set_defaults(func=_run_segment)

    return parser


def create_tiles(im, width: int, height: int, tile_size: int = 512):
    """
    Create tiles from an image.

    Args:
        im (numpy.ndarray): Input image to create tiles from.
        width (int): Width of the image.
        height (int): Height of the image.
        tile_size (int, optional): Size of the tiles. Defaults to 512.

    Returns:
        list: A list containing tiles of the input image.

    Raises:
        ValueError: If the 'tile_size' is larger than the provided image dimensions.
    """
    if tile_size > np.max(np.array(im.shape)[:-1]):
        raise ValueError("Tile size is larger than the provided image")

    tiles = []

    # Create tiles
    for x in range(0, width - tile_size, tile_size):
        for y in range(0, height - tile_size, tile_size):
            tiles.append(im[x : (x + tile_size), y : (y + tile_size)])

    return tiles


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


def cellpose_segmentation(
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
        [im],
        diameter=diameter,
        channels=channels,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
    )
    return mask


def _run_segment(args):
    """
    Wrapper for whole or tiled segmentation with cellpose.

    Args:
        args: Argument object containing input and output file paths and parameters.

    Raises:
        FileNotFoundError: If input or output directories do not exist.
    """
    logging.info("openst segmentation; running with parameters:")
    logging.info(args.__dict__)

    # Check input and output data
    check_file_exists(args.image_in)

    if not check_directory_exists(args.output_mask):
        raise FileNotFoundError("Parent directory for --output-mask does not exist")

    if args.metadata_out != "" and not check_directory_exists(args.metadata_out):
        raise FileNotFoundError("Parent directory for the metadata does not exist")

    # Set maximum image size to support large HE
    Image.MAX_IMAGE_PIXELS = args.max_image_pixels

    # Load cellpose model (path or pretrained)
    if args.model in models.MODEL_NAMES:
        model = models.Cellpose(gpu=args.gpu, model_type=args.model).cp
    else:
        check_file_exists(args.model)
        model = models.CellposeModel(gpu=args.gpu, pretrained_model=args.model)

    # Load iamge
    im = Image.open(args.image)

    # Getting width and height, we can get channels with len(im.getbands())
    width, height = im.width, im.height

    # Converting to numpy array
    im = np.array(im)

    # Prediction of mask by tile or over entire image (warning, lot of GPU memory!)
    if args.by_tiles:
        tiles = create_tiles(im, width, height, args.tile_size)
        masks = []
        for tile in tqdm(tiles):
            masks.append(
                cellpose_segmentation(
                    tile,
                    model,
                    diameter=args.diameter,
                    flow_threshold=args.flow_threshold,
                    cellprob_threshold=args.cellprob_threshold,
                )
            )
        mask_complete = assemble_from_tiles(tiles, width, height, 1, args.tile_size)
    else:
        mask_complete = cellpose_segmentation(
            tile,
            model,
            diameter=args.diameter,
            flow_threshold=args.flow_threshold,
            cellprob_threshold=args.cellprob_threshold,
        )

    dtype = np.uint8
    if mask_complete.max() >= (2**16 - 1):
        dtype = np.uint16
    elif mask_complete.max() >= (2**32 - 1):
        dtype = np.uint32
    elif mask_complete.max() >= (2**64 - 1):
        dtype = np.uint64

    if args.dilate_px > 0:
        mask_complete = expand_labels(mask_complete, distance=args.dilate_px)

    mask_complete = measure.label(mask_complete)

    Image.fromarray(mask_complete.astype(dtype)).save(args.output_mask)


if __name__ == "__main__":
    args = get_segment_parser().parse_args()
    _run_segment(args)
