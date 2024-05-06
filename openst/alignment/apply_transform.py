import logging
import json
import h5py

import numpy as np
from skimage.transform import estimate_transform as ski_estimate_transform

from openst.alignment.transformation import apply_transform
from openst.utils.file import (check_adata_structure, check_file_exists,
                               load_properties_from_adata)

def estimate_transform(model: str, src: np.ndarray, dst: np.ndarray):
    """
    TODO: write documentation
    returns sklearn transform and True/False depending on whether flip needs to be applied to the src points before tform
    """
    tform_points = ski_estimate_transform(model, src, dst)

    src_flip = np.array([src[:, 0], (src[:, 1]*-1) - ((src[:, 1]*-1).min() - (src[:, 1]).min())]).T
    tform_points_flip = ski_estimate_transform(model, src_flip, dst)

    _distance_defa = np.mean(np.linalg.norm(tform_points(src) - dst, axis=1))
    _distance_flip = np.mean(np.linalg.norm(tform_points_flip(src_flip) - dst, axis=1))

    if _distance_flip < _distance_defa:
        return tform_points_flip, True
    else:
        return tform_points, False

def apply_transform_to_coords(
    in_coords: np.ndarray,
    tile_id: np.ndarray,
    keypoints: dict,
    check_bounds: bool = False,
) -> np.ndarray:
    """
    Perform manual registration of spatial transcriptomics (STS) data with a staining image.

    Args:
        in_coords (np.ndarray): Input STS coordinates.
        total_counts (np.ndarray): Total UMI counts for each STS coordinate.
        tile_id: Identifier for each STS coordinate. During the fine registration,
                 this 'tile_id' is used to aggregate the coordinates into buckets that
                 are aligned separately. Recommended for flow-cell based STS.

    Returns:
        sts_coords_transformed (np.ndarray): Registered STS coordinates after coarse registration
    """
    sts_coords_transformed = in_coords.copy()
    if tile_id is None:
        mkpts = keypoints['all_tiles_coarse']
        mkpts_coarse0, mkpts_coarse1 = np.array(mkpts['point_src']).astype(float), np.array(mkpts['point_dst']).astype(float)
        # Preparing images and preprocessing routines
        tform_points, needs_flip = estimate_transform("similarity", mkpts_coarse1, mkpts_coarse0)
        if needs_flip:
            # this applies flipping to the coordinates
            sts_coords_transformed[..., 0] = (sts_coords_transformed[..., 0]*-1) - ((mkpts_coarse1[..., 1]*-1).min() - (mkpts_coarse1[..., 1]).min())

        sts_coords_transformed[..., :2] = apply_transform(sts_coords_transformed[..., :2], tform_points, check_bounds=check_bounds)[..., :2][..., ::-1]
    else:
        # Collect tile identifiers
        tile_codes = np.unique(tile_id.codes)

        for tile_code in tile_codes:
            if f'{tile_code}' not in keypoints.keys():
                continue
            mkpts = keypoints[f'{tile_code}']
            mkpts_fine0, mkpts_fine1 = np.array(mkpts['point_src']).astype(float), np.array(mkpts['point_dst']).astype(float)

            _t_valid_coords = tile_id.codes == tile_code
            tform_points, needs_flip = estimate_transform("similarity", mkpts_fine1, mkpts_fine0)
            if needs_flip:
                if len(mkpts_fine1) < 3:
                    logging.warn(f"Skipping tile {tile_code}: optimal transform required flipping but supported by < 3 keypoints")
                    continue
        
                # this applies flipping to the coordinates
                sts_coords_transformed[..., 0][_t_valid_coords] = ((sts_coords_transformed[..., 0]*-1) - ((mkpts_fine1[..., 1]*-1).min() - (mkpts_fine1[..., 1]).min()))[_t_valid_coords]

            sts_coords_transformed[..., :2][_t_valid_coords] = apply_transform(sts_coords_transformed[..., :2][_t_valid_coords], tform_points, check_bounds=check_bounds)[..., :2][..., ::-1]

    return sts_coords_transformed

def keypoints_json_to_dict(keypoints_json):
    keypoints_by_key = {}

    for keypoint in keypoints_json:
        if keypoint['layer'] not in keypoints_by_key.keys():
            keypoints_by_key[keypoint['layer']] = {"point_src": [], "point_dst": []}

        keypoints_by_key[keypoint['layer']]["point_src"].append(keypoint["point_src_offset_rescaled"])
        keypoints_by_key[keypoint['layer']]["point_dst"].append(keypoint["point_dst_offset_rescaled"])

    return keypoints_by_key

def load_keypoints_from_json(fname: str):
    
    with open(fname) as j:
        keypoints_json = json.load(j)['points']

    return keypoints_json_to_dict(keypoints_json)

def _run_apply_transform(args):
    # Check input and output data
    check_file_exists(args.h5_in)
    check_adata_structure(args.h5_in)

    _to_load = ["obs/total_counts"]
    if args.per_tile:
        _to_load += ["obs/tile_id"]
    _to_load += [f"{args.spatial_key_in}"]

    _tile_ids = None

    logging.info(f"Loading properties from {args.h5_in}")
    sts = load_properties_from_adata(args.h5_in, properties=_to_load)

    logging.info(f"Loading manually selected keypoints from {args.keypoints_in}")
    keypoints = load_keypoints_from_json(args.keypoints_in)
    
    logging.info(f"Applying coordinate transformation")
    # transpose spatial locations to XY coordinates
    _coords = sts[f"{args.spatial_key_in}"][:][..., ::-1]

    if args.per_tile:
        _tile_ids = sts["obs/tile_id"]
    
    _coords_transformed = apply_transform_to_coords(
        _coords,
        _tile_ids,
        keypoints
    )

    logging.info(f"Modifying coordinates in {args.h5_in}")
    with h5py.File(args.h5_in, 'r+') as adata:
        if f"{args.spatial_key_out}" in adata:
            # reconvert to YX (same axes as the images)
            adata[f"{args.spatial_key_out}"][...] = _coords_transformed[:][..., ::-1]
        else:
            adata[f"{args.spatial_key_out}"] = _coords_transformed[:][..., ::-1]

    logging.info(f"Output {args.h5_in} file was written. Finished!")

if __name__ == "__main__":
    from openst.cli import get_apply_transform_parser
    args = get_apply_transform_parser().parse_args()
    _run_apply_transform(args)
