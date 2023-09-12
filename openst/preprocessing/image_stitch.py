import argparse
import logging
import os
import subprocess

from openst.utils.file import check_directory_exists, check_file_exists


def get_image_stitch_parser():
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="stitching various image tiles into a single image (wrapper for ImageJ)",
    )

    parser.add_argument(
        "--input-dir",
        type=str,
        help="path to collection of images",
        required=True,
    )

    parser.add_argument(
        "--output-image",
        type=str,
        help="path where to save the image (must be a filename)",
        required=True,
    )

    parser.add_argument(
        "--imagej-bin",
        type=str,
        help="path to the ImageJ/Fiji executable. Must have the Grid Collection plugin available!",
        required=True,
    )

    parser.add_argument(
        "--microscope",
        type=str,
        help="microscope model or imaging strategy that was used for imaging",
        choices=["keyence"],
        required=True,
    )

    parser.add_argument(
        "--no-run",
        default=False,
        action="store_true",
        help="When specified, do not run ImageJ, but return the command line",
    )

    parser.add_argument(
        "--metadata-out",
        type=str,
        default="",
        help="Path where the metadata will be stored. If not specified, metadata is not saved.",
    )

    return parser


def setup_image_stitch_parser(parent_parser):
    """setup_image_stitch_parser"""
    parser = parent_parser.add_parser(
        "spatial_stitch",
        help="stitching STS tiles into a global coordinate system",
        parents=[get_image_stitch_parser()],
    )
    parser.set_defaults(func=_run_image_stitch)

    return parser


def _image_stitch_imagej_keyence(
    imagej_bin: str,
    input_dir: str,
    output_image: str,
):
    import binascii
    import zipfile

    macro_keyence_path = "imagej_macros/keyence_stitch.ijm"

    # Checking whether the grid file exists
    check_file_exists(os.path.join(args.input_dir, "Image.bcf"))

    # Loading the bcf metadata from keyence microscope
    archive = zipfile.ZipFile(os.path.join(input_dir, "Image.bcf"), "r")
    imgdata = archive.read("GroupFileProperty/Marker/StackList")
    GRIDX = int(binascii.hexlify(imgdata)[16:18], 16)
    GRIDY = int(binascii.hexlify(imgdata)[24:26], 16)

    if os.path.exists(output_image) and os.path.getsize(output_image) > 0:
        raise FileExistsError("The file {output_image} was already generated. Will not run!")

    logging.info(f"Grid size: X = {GRIDX}, Y = {GRIDY}")
    logging.info(f"Running macro {macro_keyence_path} with ImageJ at {imagej_bin}")

    return [
        f"{imagej_bin}",
        "--headless",
        "--console",
        "-macro",
        f"{macro_keyence_path}",
        f"{GRIDX};{GRIDY};{input_dir};{output_image}",
    ]


SUPPORTED_MICROSCOPES = {"keyence": _image_stitch_imagej_keyence}


def image_stitch_imagej(imagej_bin: str, microscope: str, input_dir: str, output_image: str):
    check_file_exists(imagej_bin)

    if microscope not in SUPPORTED_MICROSCOPES.keys():
        raise NotImplementedError(
            f"Microscope {microscope} is not supported. Only {SUPPORTED_MICROSCOPES} are supported"
        )

    return SUPPORTED_MICROSCOPES[microscope](imagej_bin, input_dir, output_image)


def _run_image_stitch(args):
    """_run_image_stitch."""
    logging.info("openst image stitching; running with parameters:")
    logging.info(args.__dict__)

    # Check input and output data
    check_file_exists(args.input_dir)

    if not check_directory_exists(args.output_image):
        raise FileNotFoundError("Parent directory for --output-mask does not exist")

    if args.metadata_out != "" and not check_directory_exists(args.metadata_out):
        raise FileNotFoundError("Parent directory for the metadata does not exist")

    cmd = image_stitch_imagej(
        imagej_bin=args.imagej_bin,
        input_dir=args.input_dir,
        output_image=args.output_image,
    )

    if args.no_run:
        logging.info("The ImageJ macro can be run with the following command")
        print(" ".join(cmd))
    else:
        _ = subprocess.run(cmd)


if __name__ == "__main__":
    args = get_image_stitch_parser().parse_args()
    _run_image_stitch(args)
