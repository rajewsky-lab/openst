import cv2
import numpy as np
from skimage.transform import resize
from skimage.filters import gaussian
import logging

def create_paired_pseudoimage(
    coords: np.ndarray,
    scale: float,
    target_size: tuple = None,
    valid_locations: np.ndarray = None,
    recenter=True,
    rescale=True,
    values=None,
    resize_method: str = 'scikit-image',
) -> dict:
    """
    Create a pseudoimage representation from input coordinates (two-dimensional), paired
    to a staining image (will set the limits and dimensions accordingly)

    Args:
        coords (np.ndarray): Input coordinates to create the pseudoimage from.
        scale (float): Scaling factor applied to rescaled coordinates.
        target_size (tuple, optional): Target size for the pseudoimage. If provided, preserves aspect ratio.
        valid_locations (np.ndarray, optional): Boolean mask indicating valid locations for creating the pseudoimage.
        recenter (bool, optional): If True, the minimum (x, y) coordinate will be offset to a new (0, 0).
        rescale (bool, optional): If True, a new scaling will be applied according to the argument 'scale'.
        values (np.ndarray, optional): When not None, will be used to populate the image intensity values.

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
    offset_factor = None

    if recenter:
        offset_factor = coords.min(axis=0)
        coords -= offset_factor

    # Rescale the coordinates to have approximately PSEUDOIMG_SIZE
    if rescale:
        rescale_factor = coords.max(axis=0).min()
        coords_rescaled = coords / rescale_factor
        coords_rescaled *= scale
        coords_rescaled_int = coords_rescaled.astype(int)
    else:
        rescale_factor = None
        coords_rescaled = coords
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

    if values is None:
        values = 1
    
    # new version to generate the pseudoimage
    # np.add.at(_sts_pseudoimage, (coords_rescaled_int[:, 0], coords_rescaled_int[:, 1]), values)
    _sts_pseudoimage, _, _ = np.histogram2d(coords_rescaled_int[:, 0], coords_rescaled_int[:, 1],
                                   bins=(dim_1, dim_2),
                                   range=[[0, dim_1], [0, dim_2]])
    #_sts_pseudoimage[coords_rescaled_int[:, 0], coords_rescaled_int[:, 1]] += values

    if resize_method == 'scikit-image':
        sts_pseudoimage = resize(_sts_pseudoimage, target_size[:2], anti_aliasing=True)
    elif resize_method == 'cv2':
        _ker = 2
        if np.array(target_size).max() > 5000:
            _ker = 6
        elif np.array(target_size).max() > 1000:
            _ker = 4
        _sts_pseudoimage = cv2.blur(_sts_pseudoimage, (_ker, _ker))
        sts_pseudoimage = cv2.resize(_sts_pseudoimage, target_size[:2][::-1], interpolation=cv2.INTER_NEAREST)
    else:
        raise NotImplementedError(f"The resize method '{resize_method}' was not implemented")
    
    sts_pseudoimage = sts_pseudoimage - sts_pseudoimage.min()
    sts_pseudoimage = ((sts_pseudoimage / sts_pseudoimage.max()) * 255).astype(np.uint8)

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


# pseudoimage for density segmentation
def recenter_points(points):
    points_roi = points.copy()
    points_roi[:, 0] = points_roi[:, 0] - points[:, 0].min()
    points_roi[:, 1] = points_roi[:, 1] - points[:, 1].min()
    return points_roi

def show_expression_on_image(points_roi,
                             render_scale: int = 1, 
                             render_sigma: float = 1.5,
                             output_resolution: float = 1):
    im_shape = points_roi.max(axis=0)

    gene_im, _, _ = np.histogram2d(points_roi[:, 0], points_roi[:, 1],
                                   bins=tuple((im_shape * (1/render_scale)).astype(int)))
    gene_im = gaussian(gene_im, render_sigma)
    gene_im = resize(gene_im, tuple((im_shape * output_resolution).astype(int)))
    return gene_im

def create_unpaired_pseudoimage(
    adata,
    spatial_coord_key: str = "obsm/spatial",
    input_resolution: float = 1,
    render_scale: float = 1,
    render_sigma: float = 1.5,
    output_resolution: float = 1,
    write_rescaled: bool = True
):
    """
    Create pseudoimage for segmentation based on RNA density (experimental feature), i.e.,
    not paired to a staining image (no cropping, no rescaling)

    Args:
        adata (ad.AnnData): Input AnnData object containing the spot-by-gene matrix.
        lims (tuple): the 2D limits for cropping the spatial coordinates, as (x_min, x_max, y_min, y_max).
        shape (tuple): the final shape of the image that will be segmented. Should lead to 1:1 aspect ratio.
        render_scale (float): rescale the coordinates by render_scale for calculating the bin image
        render_sigma (float): apply gaussian smoothing with render_sigma to the rendered pseudoimage.

    Returns:
        numpy.ndarray: pseudoimage
        numpy.ndarray: transformed points
    """

    _spatial_coords = adata[spatial_coord_key][:]
    _total_counts = adata["obs/total_counts"][:].astype(int)
    try:
        _pct_mt_counts = np.nan_to_num(adata["obs/pct_counts_mt"][:])/100
        _total_counts = (_total_counts - _total_counts * _pct_mt_counts).astype(int)
    except:
        logging.info("'pct_counts_mt' was not found; pseudoimage may contain mitochondrial counts")

    marker_filtered = recenter_points(_spatial_coords) * input_resolution
    marker_filtered_repeat = marker_filtered[np.repeat(np.arange(len(marker_filtered)), _total_counts)]

    pim = show_expression_on_image(marker_filtered_repeat, render_scale, render_sigma, output_resolution)
    logging.info(f"Created pseudoimage with {pim.shape} pixels")

    marker_filtered_scaled = marker_filtered * output_resolution

    # we need to write the transformed coordinates so they can be applied to the pseudo image
    if write_rescaled:
        _out_spatial_coord_key = f"{spatial_coord_key}_pseudoimage_scale_{render_scale}_sigma_{render_sigma}"
        
        if _out_spatial_coord_key in adata:
            adata[_out_spatial_coord_key][...] = marker_filtered_scaled
        else:
            adata[_out_spatial_coord_key] = marker_filtered_scaled

        logging.info(f"Added transformed coordinates as {_out_spatial_coord_key} to the AnnData")

    return pim, marker_filtered_scaled

def _run_pseudoimage_visualizer(args):
    import h5py
    try:
        import napari
    except ImportError:
        raise ImportError(
            "Please install napari: `pip install napari`."
        )

    from openst.utils.file import check_file_exists
    
    check_file_exists(args.h5_in)
    adata = h5py.File(args.h5_in, 'r+')
    
    im, pts = create_unpaired_pseudoimage(adata, 
                                     args.spatial_coord_key,
                                     args.input_resolution,
                                     args.render_scale,
                                     args.render_sigma,
                                     args.output_resolution)
    
    viewer = napari.Viewer()
    viewer.add_image(data=im)
    napari.run()

if __name__ == "__main__":
    from openst.cli import get_pseudoimage_parser
    args = get_pseudoimage_parser().parse_args()
    _run_pseudoimage_visualizer()