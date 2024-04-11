import argparse
import cv2
import numpy as np
from skimage.transform import resize
from skimage.filters import gaussian
import logging

def get_preview_parser():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="preview pseudoimage and images of Open-ST data",
        allow_abbrev=False,
        add_help=False,
    )
    # Input
    parser.add_argument(
        "--adata",
        type=str,
        required=True,
        help="Necessary to create the pseudoimage",
    )

    # RNA density-based pseudoimage
    parser.add_argument(
        "--spatial-coord-keys",
        type=str,
        nargs="+",
        default=None,
        help="""Path to the spatial coordinates inside the AnnData object (e.g., 'obsm/spatial').
                Can be one or many (separated by space)"""
    )

    # Staining image
    parser.add_argument(
        "--image-keys",
        type=str,
        nargs="+",
        default=None,
        help="""Path to the image to be visualized.
              Can be one or many (separated by space)"""
    )

    # Resampling before previsualizing
    parser.add_argument(
        "--spatial-coord-resampling",
        type=int,
        nargs="+",
        default=[1],
        help="""Will load every n-th point. Can be one (same for all spatial-coords)
                or many (1-to-1 mapping to the spatial-coord list)"""
    )
    parser.add_argument(
        "--image-resampling",
        type=int,
        nargs="+",
        default=[1],
        help="""Will load every n-th pixel. Can be one (same for all images)
                or many (1-to-1 mapping to the image list)"""
    )

    return parser

def setup_preview_parser(parent_parser):
    """setup_preview_parser"""
    parser = parent_parser.add_parser(
        "preview",
        help="prepreview the dataset with napari",
        parents=[get_preview_parser()],
    )
    parser.set_defaults(func=_run_preview)

    return parser

def _run_preview(args):
    import h5py
    try:
        import napari
    except ImportError:
        raise ImportError(
            "Please install napari: `pip install napari`."
        )

    from openst.utils.file import check_file_exists
    
    check_file_exists(args.adata)
    adata = h5py.File(args.adata, 'r+')

    viewer = napari.Viewer()

    if args.image_keys is not None:
        for i, _image_key in enumerate(args.image_keys):
            if len(args.image_resampling) != len(args.image_keys):
                image_resampling = args.image_resampling[0]
            else:
                image_resampling = args.image_resampling[i]

            viewer.add_image(data=adata[_image_key][::image_resampling, ::image_resampling])

    if args.spatial_coord_keys is not None:
        for i, _spatial_coord_key in enumerate(args.spatial_coord_keys):
            if len(args.spatial_coord_resampling) != len(args.spatial_coord_keys):
                spatial_coord_resampling = args.spatial_coord_resampling[0]
            else:
                spatial_coord_resampling = args.spatial_coord_resampling[i]

            viewer.add_points(data=adata[_spatial_coord_key][::spatial_coord_resampling])
    
    napari.run()

if __name__ == "__main__":
    args = get_preview_parser().parse_args()
    _run_preview()