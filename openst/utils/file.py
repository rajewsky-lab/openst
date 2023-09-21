import logging
import os
import pickle
import shutil
from pathlib import Path
from typing import Union

import h5py
from anndata import AnnData
from anndata._io.specs import read_elem


def save_pickle(obj, file_path):
    """
    Save an object to a pickle file.

    Parameters:
        obj (any): The object to be saved.
        file_path (str): The path to the pickle file.

    Returns:
        None
    """
    with open(file_path, "wb") as f:
        pickle.dump(obj, f)


def load_pickle(file_path):
    """
    Load an object from a pickle file.

    Parameters:
        file_path (str): The path to the pickle file.

    Returns:
        any: The loaded object.
    """
    with open(file_path, "rb") as f:
        obj = pickle.load(f)
    return obj


def check_file_exists(f, exception=True):
    """
    Check whether the file exists.

    Args:
        f (str): Path to the input file.

    Raises:
        FileNotFoundError: If the file does not exist.
    """

    if not os.path.exists(f):
        if exception:
            raise FileNotFoundError(f"The file '{f}' does not exist.")
        else:
            return False

    return True


def check_directory_exists(path):
    """
    Check if a file exists, or if its parent directory exists.

    Parameters:
        path (str): Path to the file or directory.

    Returns:
        bool: True if the parent directory exists or if the file exists, False otherwise.
    """
    if os.path.isdir(path):
        return os.path.exists(path)
    else:
        parent_directory = os.path.dirname(path)
        # handle file created in the same directory
        if parent_directory == "":
            return True
        return os.path.exists(parent_directory)


def check_adata_structure(f):
    """
    Check the validity of the input AnnData file.

    Args:
        f (str): Path to the input AnnData file.

    Raises:
        KeyError: If required properties are not found in the file.
    """

    with h5py.File(f, "r") as file:
        if "obsm/spatial" not in file:
            raise KeyError("The AnnData file does not have the 'obsm/spatial' property.")

        if "obs/puck_id" not in file:
            raise KeyError("The AnnData file does not have the 'obs/puck_id' property.")

        if "obs/total_counts" not in file:
            raise KeyError("The AnnData file does not have the 'obs/total_counts' property.")

        if "spatial_aligned" in file:
            logging.warn("The AnnData file has a 'spatial_aligned' layer")


def load_properties_from_adata(f: Union[str, AnnData], properties: list = ["obsm/spatial"], backed: bool=False) -> dict:
    """
    Load specified properties from an AnnData file (h5py format).

    Args:
        f (str): Path to the AnnData h5py file.
        properties (list, optional): List of property paths to load from the file.
        backed (bool, optional): If True, data will not be read into memory.

    Returns:
        dict: A dictionary containing the loaded properties.
            - For each property path specified in the 'properties' list:
                * The dictionary key is the property path.
                * The value is the corresponding parsed property data.

    Notes:
        - This function loads specified properties from an AnnData h5py file.
        - The 'properties' list should consist of property paths within the file.
        - Returns a dictionary where keys are property paths and values are the loaded data.
    """

    parsed_properties = {}

    if isinstance(f, AnnData):
        for p in properties:
            parsed_properties[p] = read_elem(f[p])
    elif isinstance(f, str):
            if backed:
                _f = h5py.File(f)
                for p in properties:
                    parsed_properties[p] = _f[p]
            else:
                with h5py.File(f) as _f:
                    for p in properties:
                        parsed_properties[p] = read_elem(_f[p])
    else:
        raise TypeError("Type of 'f' is incorrect. It needs to be an AnnData or str object.")

    return parsed_properties


def check_obs_unique(adata: AnnData, obs_key: str = "tile_id") -> bool:
    """
    Check if the values in a specified observation key in an AnnData object are unique.

    Args:
        adata (AnnData): AnnData object to check for unique observations.
        obs_key (str, optional): The name of the observation key to check for uniqueness. Defaults to "tile_id".

    Returns:
        bool: True if the specified observation key has unique values, False otherwise.

    Raises:
        ValueError: If the specified observation key exists in the AnnData object but is not unique.
    """
    return adata.obs[obs_key].nunique() == 1


def copytree2(source: str, dest: str) -> str:
    """
    Recursively copy the contents of a source directory to a destination directory.

    Args:
        source (str): The source directory to be copied.
        dest (str): The destination directory where the contents will be copied to.

    Returns:
        str: The path to the destination directory where the contents were copied.

    Notes:
        - This function creates the destination directory and its parent directories if they do not exist.
        - It checks if the source and destination directories already exist and have the same size.
          If so, it skips copying.
        - If the source and destination directories differ in size or do not exist, it performs a recursive copy.

    """
    Path(dest).mkdir(parents=True, exist_ok=True)
    dest_dir = os.path.join(dest, os.path.basename(source))
    if os.path.exists(dest_dir) and os.path.getsize(dest_dir) == os.path.getsize(source):
        print("The directory {OUTFILE} was already copied. Skipping!")
    else:
        shutil.copytree(source, dest_dir, dirs_exist_ok=True)
    return dest_dir
