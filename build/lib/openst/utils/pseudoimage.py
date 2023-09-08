import numpy as np
from skimage.transform import resize

def create_pseudoimage(
    coords: np.ndarray,
    scale: float,
    target_size: tuple = None,
    valid_locations: np.ndarray = None,
) -> dict:
    """
    Create a pseudoimage representation from input coordinates (two-dimensional).

    Args:
        coords (np.ndarray): Input coordinates to create the pseudoimage from.
        scale (float): Scaling factor applied to rescaled coordinates.
        target_size (tuple, optional): Target size for the pseudoimage. If provided, preserves aspect ratio.
        valid_locations (np.ndarray, optional): Boolean mask indicating valid locations for creating the pseudoimage.

    Returns:
        dict: A dictionary containing the pseudoimage and related metadata.
            - 'pseudoimage' (np.ndarray): The generated pseudoimage.
            - 'rescaling_factor' (float): Scaling ratio used for rescaling the pseudoimage.
            - 'rescale_factor' (float): Scaling factor applied to the input coordinates.
            - 'offset_factor' (np.ndarray): Offset applied to the input coordinates.
            - 'target_size' (tuple): Target size for the pseudoimage (if provided).
            - 'scale' (float): Scaling factor applied to the rescaled coordinates.
            - 'valid_locations' (np.ndarray): Boolean mask indicating valid locations for creating the pseudoimage.
            - 'coords_rescaled' (np.ndarray): Rescaled coordinates used for generating the pseudoimage.

    Notes:
        - The function generates a pseudoimage by rescaling the input coordinates and populating a grid.
        - The rescaled coordinates are used to calculate the scaling ratio and apply transformations.
        - If target_size is provided, the pseudoimage is resized to match the target size.
    """

    if type(coords) is not np.ndarray:
        raise TypeError("'coords' is expected to be of type np.ndarray")
    elif coords.ndim == 2 and coords.shape[1] != 2:
        raise ValueError("Input coordinates must be two-dimensional")
    elif coords.ndim != 2:
        raise ValueError("'coords' array does not have the expected number of dimensions")

    # Apply coordinate transformation (zero-rescaling)
    coords = coords.copy()
    offset_factor = coords.min(axis=0)
    coords -= offset_factor

    # Rescale the coordinates to have approximately PSEUDOIMG_SIZE
    rescale_factor = coords.max(axis=0).min()
    coords_rescaled = coords / rescale_factor
    coords_rescaled *= scale
    coords_rescaled_int = coords_rescaled.astype(int)

    dim_1, dim_2 = coords_rescaled_int.max(axis=0)

    # Preserve aspect ratio given a target_size
    if target_size is not None:
        dim_2_prop = dim_1 / (target_size[0] / target_size[1])
        dim_1_prop = dim_2 / (target_size[1] / target_size[0])

        if dim_2_prop < dim_2:
            dim_1 = int(dim_1_prop) + 1
            dim_2 = int(dim_1_prop / (target_size[0] / target_size[1])) + 1
        elif dim_1_prop < dim_1:
            dim_2 = int(dim_2_prop) + 1
            dim_1 = int(dim_2_prop / (target_size[1] / target_size[0])) + 1

    _sts_pseudoimage = np.zeros((dim_1 + 1, dim_2 + 1))

    if valid_locations is not None:
        coords_rescaled_int = coords_rescaled_int[valid_locations]

    _sts_pseudoimage[coords_rescaled_int[:, 0], coords_rescaled_int[:, 1]] += 1

    sts_pseudoimage = resize(_sts_pseudoimage, target_size[:2], anti_aliasing=True)

    # calculate scaling ratio from the initially transformed points;
    # use these points for applying the transform matrix (not integer)
    rescaling_factor = sts_pseudoimage.shape[0] / _sts_pseudoimage.shape[0]

    pseudoimage_and_metadata = {
        "pseudoimage": sts_pseudoimage,
        "rescaling_factor": rescaling_factor,
        "rescale_factor": rescale_factor,
        "offset_factor": offset_factor,
        "target_size": target_size,
        "scale": scale,
        "valid_locations": valid_locations,
        "coords_rescaled": coords_rescaled,
    }

    return pseudoimage_and_metadata