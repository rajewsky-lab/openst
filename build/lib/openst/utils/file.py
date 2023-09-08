import h5py
import logging
import os
import pickle

from anndata import read_h5ad
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


def check_file_exists(f):
    """
    Check whether the file exists.

    Args:
        f (str): Path to the input file.

    Raises:
        FileNotFoundError: If the file does not exist.
    """

    if not os.path.exists(f):
        raise FileNotFoundError(f"The file '{f}' does not exist.")
    

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


def load_properties_from_adata(f: str, properties: list = ["obsm/spatial"]) -> dict:
    """
    Load specified properties from an AnnData file (h5py format).

    Args:
        f (str): Path to the AnnData h5py file.
        properties (list, optional): List of property paths to load from the file.

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

    with h5py.File(f) as f:
        for p in properties:
            parsed_properties[p] = read_elem(f[p])

    return parsed_properties


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
        return os.path.exists(parent_directory)