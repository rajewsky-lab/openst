import logging
import os
import pickle
import shutil
from pathlib import Path


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


def check_file_exists(f, exception=True) -> bool:
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


def check_directory_exists(path, exception=False) -> bool:
    """
    Check if a file exists, or if its parent directory exists.

    Parameters:
        path (str): Path to the file or directory.

    Returns:
        bool: True if the parent directory exists or if the file exists, False otherwise.
    """
    _ret_val = False
    if os.path.isdir(path):
        _ret_val = os.path.exists(path)
    else:
        path = os.path.dirname(path)
        # handle file created in the same directory
        if path == "":
            _ret_val = True
        else:
            _ret_val = os.path.exists(path)
    
    if exception and not _ret_val:
        raise FileNotFoundError(f"The directory '{path}' does not exist")
    
    return _ret_val


def check_adata_structure(f):
    """
    Check the validity of the input Open-ST h5 object.

    Args:
        f (str): Path to the input Open-ST h5 object.

    Raises:
        KeyError: If required properties are not found in the file.
    """
    import h5py

    with h5py.File(f, "r") as file:
        if "obsm/spatial" not in file:
            raise KeyError("The Open-ST h5 object does not have the 'obsm/spatial' property.")

        if "obs/tile_id" not in file:
            raise KeyError("The Open-ST h5 object does not have the 'obs/tile_id' property.")

        if "obs/total_counts" not in file:
            raise KeyError("The Open-ST h5 object does not have the 'obs/total_counts' property.")

        if "spatial_aligned" in file:
            logging.warn("The Open-ST h5 object has a 'spatial_aligned' layer")


def load_properties_from_adata(f, properties: list = ["obsm/spatial"], backed: bool=False) -> dict:
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
    import h5py
    from anndata import AnnData
    from anndata._io.specs import read_elem

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


def check_obs_unique(adata, obs_key: str = "tile_id") -> bool:
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

def get_package_path() -> str:
    """Get the absolute path of the directory containing the current Python package.

    Returns:
        str: Absolute path of the directory containing the current Python package.
    """
    import openst
    return os.path.dirname(os.path.abspath(openst.__file__))

def get_absolute_package_path(relative_path) -> str:
    """
    Get the absolute path by concatenating the package path and the relative path.

    Args:
        relative_path (str): Relative path from the package directory.

    Returns:
        str: Absolute path.
    """
    package_path = get_package_path()
    return os.path.join(package_path, relative_path)

def h5_to_dict(adata) -> dict:
    """
    Recursively converts an h5py.Group object and its nested datasets into a nested dictionary structure.

    Parameters:
        adata (h5py.Group): An h5py Group object to be converted.

    Returns:
        dict: A nested dictionary representing the structure of the h5py Group object.
            Leaf nodes contain strings representing the type and shape (if applicable) of the datasets.
            Non-leaf nodes contain nested dictionaries representing their child groups and datasets.

    Notes:
        - Leaf nodes in the resulting dictionary contain strings formatted as "{type}_{shape}".
          If the dataset has no shape attribute (e.g., scalar dataset), shape will be None.
          Example: "<class 'numpy.ndarray'>_(10,)"
        - Non-leaf nodes in the resulting dictionary contain nested dictionaries
          representing their child groups and datasets.
    """
    import h5py

    result = {}
    for key, value in adata.items():
        if isinstance(value, h5py.Group):
            result[key] = h5_to_dict(value)
        else:
            dataset_type = str(type(value))
            dataset_shape = value.shape if hasattr(value, 'shape') else None
            result[key] = f"{dataset_type}_{dataset_shape}"
    return result

def write_key_to_h5(adata, key, data, delete_before=False):
    if key in adata and not delete_before:
        adata[key][:] = data
    elif key in adata and delete_before:
        del adata[key]
    else:
        adata[key] = data

def download_url_to_file(url, dst, progress=True):
    r"""Download object at the given URL to a local path.
            Thanks to torch & cellpose
    Args:
        url (string): URL of the object to download
        dst (string): Full path where object will be saved, e.g. `/tmp/temporary_file`
        progress (bool, optional): whether or not to display a progress bar to stderr
            Default: True
    """

    from urllib.request import urlopen
    import tempfile
    from tqdm import tqdm

    file_size = None
    import ssl
    ssl._create_default_https_context = ssl._create_unverified_context
    u = urlopen(url)
    meta = u.info()
    if hasattr(meta, "getheaders"):
        content_length = meta.getheaders("Content-Length")
    else:
        content_length = meta.get_all("Content-Length")
    if content_length is not None and len(content_length) > 0:
        file_size = int(content_length[0])
    # We deliberately save it in a temp file and move it after
    dst = os.path.expanduser(dst)
    dst_dir = os.path.dirname(dst)
    f = tempfile.NamedTemporaryFile(delete=False, dir=dst_dir)
    try:
        with tqdm(total=file_size, disable=not progress, unit="B", unit_scale=True,
                  unit_divisor=1024) as pbar:
            while True:
                buffer = u.read(8192)
                if len(buffer) == 0:
                    break
                f.write(buffer)
                pbar.update(len(buffer))
        f.close()
        shutil.move(f.name, dst)
    finally:
        f.close()
        if os.path.exists(f.name):
            os.remove(f.name)
