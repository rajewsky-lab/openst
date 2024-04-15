import logging
from typing import List, Union
from tqdm import tqdm

import numpy as np
from anndata import AnnData, concat, read_h5ad

from openst.utils.file import (check_directory_exists, check_file_exists,
                               check_obs_unique)

DEFAULT_REGEX_TILE_ID = "(L[1-4][a-b]_tile_[1-2][0-7][0-9][0-9])"

def _transform_tile(tile: AnnData, tiles_transform: dict):
    """
    Transform the spatial coordinates of an AnnData object based on a transformation dictionary.

    Args:
        tile (AnnData): AnnData object containing spatial coordinates to be transformed.
        tiles_transform (dict): A dictionary containing transformation information.

    Returns:
        AnnData: Transformed AnnData object with updated spatial coordinates.

    Raises:
        ValueError: If the 'tile_id' exists in the AnnData object but is not unique.
    """
    if not check_obs_unique(tile, "tile_id"):
        raise ValueError("tile_id exists in AnnData object but are not unique")

    _tile_id = np.unique(tile.obs["tile_id"])[0]
    x_ofs = tiles_transform["x_offset"][_tile_id]
    y_ofs = tiles_transform["y_offset"][_tile_id]
    tile.obsm["spatial"][:, 0] += x_ofs
    tile.obsm["spatial"][:, 1] += y_ofs

    return tile


def create_spatial_stitch(
    tile: AnnData,
    tile_transform: dict,
    reset_index: bool = True,
    transform: bool = True,
):
    """
    Create a tile collection from an AnnData object, optionally resetting index and applying spatial transformation.

    Args:
        tile (AnnData): AnnData object to create a tile collection from.
        tile_transform (dict): A dictionary containing transformation information.
        reset_index (bool, optional): Reset the index of the AnnData object. Defaults to True.
        transform (bool, optional): Apply spatial transformation. Defaults to True.

    Returns:
        AnnData: Tile collection AnnData object.

    Notes:
        - This function can reset the index and apply spatial transformation to create a tile collection.
    """
    if reset_index:
        tile.obs_names = tile.obs_names.astype(str) + ":" + tile.obs["tile_id"].astype(str)

    if transform:
        tile = _transform_tile(tile, tile_transform)

    return tile


def parse_tile_id_from_path(f: str, tile_id_regex: str = DEFAULT_REGEX_TILE_ID):
    """
    Parse the tile ID from a file path using a regular expression.

    Args:
        f (str): File path.
        tile_id_regex (str, optional): Regular expression pattern for extracting the tile ID.
                                       Defaults to DEFAULT_REGEX_TILE_ID.

    Returns:
        str: Extracted tile ID from the file path.
    """
    import os
    import re

    bname = os.path.basename(f)
    tile_id = re.findall(rf"{tile_id_regex}", bname)

    if len(tile_id) > 1:
        logging.warn("Found more than one tile_id in the path. First one (index 0) will be used.")

    tile_id = tile_id[0]

    return tile_id


def read_tiles_to_list(
    f: Union[str, List[str]],
    tile_id: Union[int, List[int], None] = None,
    tile_id_regex: str = DEFAULT_REGEX_TILE_ID,
    tile_id_key: str = "tile_id",
):
    """
    Read tile data from one or more files into a list of AnnData objects.

    Args:
        f (Union[str, List[str]]): File path or list of file paths to read tiles from.
        tile_id (Union[int, List[int], None], optional): Tile ID or list of tile IDs. Defaults to None.
        tile_id_regex (str, optional): Regular expression pattern for extracting tile IDs from file paths.
                                       Defaults to DEFAULT_REGEX_TILE_ID.
        tile_id_key (str, optional): Observation key name for tile IDs in the AnnData object. Defaults to "tile_id".

    Returns:
        List[AnnData]: List of AnnData objects representing tiles.
    """
    if type(f) is str:
        f = [f]

    if tile_id is not None and type(tile_id) is str:
        tile_id = [tile_id]
    elif type(tile_id) is list and len(tile_id) != len(f):
        raise ValueError(f"""You must provide {len(tile_id)} items in --tile-id,
                           one per file in --tiles (currently provides {len(f)})""")

    tiles = []

    for i, f in tqdm(enumerate(f)):
        _f_obj = read_h5ad(f)

        if "spatial" not in _f_obj.obsm.keys():
            raise ValueError(f"Could not find valid .obsm['spatial'] data in {f}")
        if tile_id_key not in _f_obj.obs.keys() and tile_id is None:
            if tile_id is None:
                _tile_id = parse_tile_id_from_path(f, tile_id_regex=tile_id_regex)
                if len(_tile_id) == 0:
                    raise ValueError(
                        f"Could not find a 'tile_id' for tile {f} with the regular expression {tile_id_regex}"
                    )
            else:
                _tile_id = tile_id[i]

            _f_obj.obs[tile_id_key] = _tile_id
        elif tile_id is not None:
            _tile_id = tile_id[i]
            _f_obj.obs[tile_id_key] = _tile_id

        if tile_id_key != "tile_id":
            _f_obj.obs["tile_id"] = _tile_id

        if not check_obs_unique(_f_obj, "tile_id"):
            raise ValueError(f"'tile_id' exists in Open-ST h5 object but contains more than one unique value in tile {f}")

        tiles.append(_f_obj)

    return tiles


def parse_tile_coordinate_system_file(f: str):
    """
    Parse a tile coordinate system file into a dictionary.

    Args:
        f (str): File path to the tile coordinate system file.

    Returns:
        dict: Dictionary containing tile coordinate system information.
    """
    import pandas as pd

    cs = pd.read_csv(f, sep="[,|\t]", engine="python")

    _tile_id_col = cs.columns[0]
    cs = cs.set_index(_tile_id_col)

    cs = cs.loc[~cs.index.duplicated(keep="first")]

    return cs.to_dict(orient="dict")


def merge_tiles_to_collection(
    tiles: List[str],
    tile_id: List[str],
    tile_coordinates: str,
    tile_id_regex: str = None,
    tile_id_key: str = "tile_id",
    no_reset_index: bool = False,
    no_transform: bool = False,
    merge_output: str = "same",
    join_output: str = "inner",
):
    """
    Merge multiple tiles into a single tile collection AnnData object.

    Args:
        tiles (List[str]): List of file paths containing tile data.
        tile_id (List[str]): List of tile IDs.
        tile_coordinates (str): File path to the tile coordinate system file.
        tile_id_regex (str, optional): Regular expression pattern for extracting tile IDs from file paths.
                                       Defaults to None.
        tile_id_key (str, optional): Observation key name for tile IDs in the AnnData object. Defaults to "tile_id".
        no_reset_index (bool, optional): Skip resetting the index of the AnnData objects. Defaults to False.
        no_transform (bool, optional): Skip applying spatial transformation. Defaults to False.
        merge_output (str, optional): Merge method for AnnData objects. Defaults to "same".
        join_output (str, optional): Join method for AnnData objects. Defaults to "inner".

    Returns:
        AnnData: Merged tile collection AnnData object.
    """
    tile_transform = parse_tile_coordinate_system_file(tile_coordinates)

    tiles_list = read_tiles_to_list(tiles, tile_id, tile_id_regex, tile_id_key)

    spatial_stitch_list = []

    for tile in tiles_list:
        spatial_stitch_list += [create_spatial_stitch(tile, tile_transform, ~no_reset_index, ~no_transform)]

    spatial_stitch = concat(spatial_stitch_list, merge=merge_output, join=join_output)

    spatial_stitch.uns = {np.unique(tile.obs[tile_id_key])[0]: tile.uns for tile in spatial_stitch_list}

    return spatial_stitch


def _run_spatial_stitch(args):
    """_run_spatial_stitch."""
    # Check input and output data
    if type(args.tiles) is str:
        check_file_exists(args.tiles)
    else:
        for t in args.tiles:
            check_file_exists(t)

    if not check_directory_exists(args.h5_out):
        raise FileNotFoundError("Parent directory for --h5-out does not exist")

    if args.metadata != "" and not check_directory_exists(args.metadata):
        raise FileNotFoundError("Parent directory for the metadata does not exist")

    logging.info(f"Loading {len(args.tiles)} tiles, correcting coordinates and merging")
    spatial_stitch = merge_tiles_to_collection(
        tiles=args.tiles,
        tile_id=args.tile_id,
        tile_coordinates=args.tile_coordinates,
        tile_id_regex=args.tile_id_regex,
        tile_id_key=args.tile_id_key,
        no_reset_index=args.no_reset_index,
        no_transform=args.no_transform,
        merge_output=args.merge_output,
        join_output=args.join_output,
    )

    logging.info(f"Writing merged tiles into {args.h5_out}")
    spatial_stitch.write_h5ad(args.h5_out)


if __name__ == "__main__":
    from openst.cli import get_spatial_stitch_parser
    args = get_spatial_stitch_parser().parse_args()
    _run_spatial_stitch(args)
