import argparse
import logging
import json
import h5py
import shutil

import numpy as np
from anndata import read_h5ad
from skimage.transform import estimate_transform

from openst.alignment.transformation import apply_transform
from openst.utils.file import (check_adata_structure, check_directory_exists,
                               check_file_exists, load_properties_from_adata)

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
        "--keypoints-json",
        type=str,
        required=True,
        help="Path to the json file containing keypoints",
    )
    parser.add_argument(
        "--h5-in",
        type=str,
        required=True,
        help="Path to the input h5ad file containing spatial coordinates",
    )
    parser.add_argument(
        "--h5-out",
        type=str,
        default="",
        help="""Path where the h5ad file will be saved after alignment.
        If not indicated, data is written in place at --h5-in""",
    )
    parser.add_argument(
        "--per-tile",
        action="store_true",
        help="If selected, individual transformations per tile are estimated from they keypoints",
    )
    parser.add_argument(
        "--spatial-key-in",
        type=str,
        default="spatial_pairwise_aligned_coarse",
        help="""The name of the `obsm` variable where the transformed coordinates will be read from""",
    )
    parser.add_argument(
        "--spatial-key-out",
        type=str,
        default="spatial_pairwise_aligned_fine",
        help="""The name of the `obsm` variable where the transformed coordinates will be written""",
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

def apply_transform_to_coords(
    in_coords: np.ndarray,
    tile_id: np.ndarray,
    keypoints: dict,
) -> (np.ndarray, np.ndarray):
    """
    Perform manual registration of spatial transcriptomics (STS) data with a staining image.

    Args:
        in_coords (np.ndarray): Input STS coordinates.
        total_counts (np.ndarray): Total UMI counts for each STS coordinate.
        tile_id: Identifier for each STS coordinate. During the fine registration,
                 this 'tile_id' is used to aggregate the coordinates into buckets that
                 are aligned separately. Recommended for flow-cell based STS.

    Returns:
        tuple: A tuple containing four elements:
            - sts_coords_transformed (np.ndarray): Registered STS coordinates after coarse registration
    """
    sts_coords_transformed = in_coords.copy()
    if tile_id is None:
        mkpts = keypoints['all_tiles_coarse']
        mkpts_coarse0, mkpts_coarse1 = np.array(mkpts['point_src']).astype(float), np.array(mkpts['point_dst']).astype(float)
        # Preparing images and preprocessing routines
        tform_points = estimate_transform("similarity", mkpts_coarse1, mkpts_coarse0)
        sts_coords_transformed[..., :2] = apply_transform(in_coords[..., :2], tform_points, check_bounds=False)[..., :2][..., ::-1]
    else:
        # Collect tile identifiers
        tile_codes = np.unique(tile_id.codes)

        for tile_code in tile_codes:
            if f'{tile_code}' not in keypoints.keys():
                continue
            mkpts = keypoints[f'{tile_code}']
            mkpts_fine0, mkpts_fine1 = np.array(mkpts['point_src']).astype(float), np.array(mkpts['point_dst']).astype(float)

            _t_valid_coords = tile_id.codes == tile_code
            tform_points = estimate_transform("similarity", mkpts_fine1, mkpts_fine0)
            sts_coords_transformed[..., :2][_t_valid_coords] = apply_transform(in_coords[..., :2][_t_valid_coords], tform_points, check_bounds=False)[..., :2][..., ::-1]

    return sts_coords_transformed

def load_keypoints_from_json(fname: str):
    keypoints_by_key = {}
    with open(fname) as j:
        _keypoints_dict = json.load(j)['points']

    for keypoint in _keypoints_dict:
        if keypoint['layer'] not in keypoints_by_key.keys():
            keypoints_by_key[keypoint['layer']] = {"point_src": [], "point_dst": []}

        keypoints_by_key[keypoint['layer']]["point_src"].append(keypoint["point_src_offset_rescaled"])
        keypoints_by_key[keypoint['layer']]["point_dst"].append(keypoint["point_dst_offset_rescaled"])

    return keypoints_by_key

def _run_manual_pairwise_aligner(args):
    logging.info("open-ST pairwise alignment; running with parameters:")
    logging.info(args.__dict__)

    # Check input and output data
    check_file_exists(args.h5_in)
    check_adata_structure(args.h5_in)

    if args.h5_out == "" and not check_directory_exists(args.h5_out):
        raise FileNotFoundError("Parent directory for --h5-out does not exist")
    
    if args.h5_out == "":
        args.h5_out = args.h5_in

    _to_load = ["obs/total_counts"]
    if args.per_tile:
        _to_load += ["obs/tile_id"]
    _to_load += [f"obsm/{args.spatial_key_in}"]

    _tile_ids = None

    logging.info(f"Loading properties from {args.h5_in}")
    sts = load_properties_from_adata(args.h5_in, properties=_to_load)

    logging.info(f"Loading manually selected keypoints from {args.keypoints_json}")
    keypoints = load_keypoints_from_json(args.keypoints_json)
    
    logging.info(f"Applying coordinate transformation")
    # transpose spatial locations to XY coordinates
    _coords = sts[f"obsm/{args.spatial_key_in}"][:][..., ::-1]

    if args.per_tile:
        _tile_ids = sts["obs/tile_id"]
    
    _coords_transformed = apply_transform_to_coords(
        _coords,
        _tile_ids,
        keypoints
    )

    if check_file_exists(args.h5_out, exception = False):
        logging.info(f"The output file exists at {args.h5_out}")  
    else:
        logging.info(f"Creating new {args.h5_out} from {args.h5_in}")  
        shutil.copy(args.h5_in, args.h5_out)

    logging.info(f"Modifying coordinates in {args.h5_out}")
    with h5py.File(args.h5_out, 'r+') as adata:
        if f"obsm/{args.spatial_key_out}" in adata:
            # reconvert to YX (same axes as the images)
            adata[f"obsm/{args.spatial_key_out}"][...] = _coords_transformed[:][..., ::-1]
        else:
            adata[f"obsm/{args.spatial_key_out}"] = _coords_transformed[:][..., ::-1]

    logging.info(f"Output {args.h5_out} file was written. Finished!")

if __name__ == "__main__":
    args = get_manual_pairwise_aligner_parser().parse_args()
    _run_manual_pairwise_aligner(args)
