import dask
import numpy as np
import logging

from dask_image.ndmorph import binary_dilation as dask_binary_dilation
from dask_image.ndfilters import gaussian as dask_gaussian
from scipy.ndimage import binary_dilation as scipy_binary_dilation
from skimage.filters import threshold_otsu
from skimage.filters import gaussian as skimage_gaussian
from skimage.color import rgb2hsv as skimage_rgb2hsv

HI_NUM_ITER = 100

def binary_fill_holes(input, structure=None, output=None, origin=0):
    """
    Fill the holes in binary objects, compatible with dask.


    Args:
        input : array_like
            N-D binary array with holes to be filled
        structure : array_like, optional
            Structuring element used in the computation; large-size elements
            make computations faster but may miss holes separated from the
            background by thin regions. The default element (with a square
            connectivity equal to one) yields the intuitive result where all
            holes in the input have been filled.
        output : ndarray, optional
            Array of the same shape as input, into which the output is placed.
            By default, a new array is created.
        origin : int, tuple of ints, optional
            Position of the structuring element.

    Returns:
        out : ndarray
            Transformation of the initial image `input` where holes have been
            filled.

    Notes:
        This was adapted from skimage

    """
    mask = np.logical_not(input)
    if type(input) is dask.array.Array:
        tmp = dask.array.zeros(mask.shape, bool)
        # TODO: dask does not support iterations until convergence, we set to a high number
        output = dask_binary_dilation(tmp, structure=structure, iterations=HI_NUM_ITER, mask=mask, border_value=1,
                                origin=origin, brute_force=False)
    else:
        tmp = np.zeros(mask.shape, bool)
        output = scipy_binary_dilation(tmp, structure, -1, mask, None, 1,
                                origin)
    
    np.logical_not(output, output)
    return output


def dask_threshold_otsu(image=None, nbins=256):
    """Return threshold value based on Otsu's method.
    WARNING: please only use with dask! Otherwise, skimage method is preferred

    Either image or hist must be provided. If hist is provided, the actual
    histogram of the image is ignored.

    Parameters:
        image : (N, M[, ..., P]) ndarray, optional
            Grayscale input image.
        nbins : int, optional
            Number of bins used to calculate histogram. This value is ignored for
            integer arrays.
        hist : array, or 2-tuple of arrays, optional
            Histogram from which to determine the threshold, and optionally a
            corresponding array of bin center intensities. If no hist provided,
            this function will compute it from the image.


    Returns:
        threshold : float
            Upper threshold value. All pixels with an intensity higher than
            this value are assumed to be foreground.

    Notes:
        Adapted from skimage; the bin calculation is different, but suffices for this purpose
        Please only use with dask!

    """
    if image is not None and image.ndim > 2 and image.shape[-1] in (3, 4):
        logging.warn(f'threshold_otsu is expected to work correctly only for '
             f'grayscale images; image shape {image.shape} looks like '
             f'that of an RGB image.')

    # Check if the image has more than one intensity value; if not, return that
    # value
    if image is not None:
        first_pixel = image.ravel()[0]
        if np.all(image == first_pixel):
            return first_pixel
        
    data_min, data_max = image.min(), image.max()
    bin_edges = np.linspace(data_min, data_max, nbins + 1) #nbins + 1

    counts, bin_edges = np.histogram(image.ravel(), bin_edges)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # class probabilities for all possible thresholds
    weight1 = np.cumsum(counts)
    weight2 = np.cumsum(counts[::-1])[::-1]
    # class means for all possible thresholds
    mean1 = np.cumsum(counts * bin_centers) / weight1
    mean2 = (np.cumsum((counts * bin_centers)[::-1]) / weight2[::-1])[::-1]

    # Clip ends to align class 1 and class 2 variables:
    # The last value of ``weight1``/``mean1`` should pair with zero values in
    # ``weight2``/``mean2``, which do not exist.
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2

    idx = np.argmax(variance12)
    threshold = bin_centers[idx]

    return threshold


def _prepare_colorarray(arr, channel_axis=-1):
    """Check the shape of the array and convert it to
    floating point representation.

    This is a simplification from the skimage implementation of the same function
    """
    if arr.shape[channel_axis] != 3:
        msg = (f'the input array must have size 3 along `channel_axis`, '
               f'got {arr.shape}')
        raise ValueError(msg)

    return arr.astype(float)


def dask_rgb2hsv(rgb, *, channel_axis=-1):
    """RGB to HSV color space conversion.

    Parameters:
        rgb : (..., 3, ...) array_like
            The image in RGB format. By default, the final dimension denotes
            channels.
        channel_axis : int, optional
            This parameter indicates which axis of the array corresponds to
            channels.

    Returns:
        out : (..., 3, ...) ndarray
            The image in HSV format. Same dimensions as input.

    Notes:
        This was adapted from skimage
    """
    input_is_one_pixel = rgb.ndim == 1
    if input_is_one_pixel:
        raise ValueError("Wrong dimensions!")

    arr = _prepare_colorarray(rgb, channel_axis=-1)
    out = np.empty_like(arr)

    # -- V channel
    out_v = arr.max(-1)

    # -- S channel
    delta = np.ptp(arr,-1)
    # Ignore warning for zero divided by zero
    out_s = delta / out_v
    out_s = np.where(delta == 0., 0., out_s)

    out.compute_chunk_sizes()
    arr.compute_chunk_sizes()
    delta.compute_chunk_sizes()

    if type(rgb) is dask.array.Array:
        out_h = np.zeros_like(out_s)
        logging.warn("Hue channel not implemented when using dask array")
    else:
        # -- H channel
        # red is max
        idx = (arr[..., 0] == out_v)
        out[idx, 0] = (arr[idx, 1] - arr[idx, 2]) / delta[idx]

        # green is max
        idx = (arr[..., 1] == out_v)
        out[idx, 0] = 2. + (arr[idx, 2] - arr[idx, 0]) / delta[idx]

        # blue is max
        idx = (arr[..., 2] == out_v)
        out[idx, 0] = 4. + (arr[idx, 0] - arr[idx, 1]) / delta[idx]
        out_h = (out[..., 0] / 6.) % 1.
        out_h[delta == 0.] = 0.

    # TODO: implement hue so it is compatible with dask

    # -- output
    out[..., 0] = out_h
    out[..., 1] = out_s
    out[..., 2] = out_v

    # # remove NaN
    out[np.isnan(out)] = 0

    if input_is_one_pixel:
        out = np.squeeze(out, axis=0)

    return out

def mask_tissue(
    image: np.ndarray,
    keep_black_background: bool = False,
    mask_gaussian_blur: float = 5,
    return_hsv: bool = True
):
    """
    Prepare an image for feature matching by applying a series of image processing steps.

    Args:
        image (np.ndarray): Input RGB image to be prepared.
        keep_black_background (bool, optional): If True, keeps background pixels as black. Only if mask_tissue=True
        mask_gaussian_blur (float, optional): Standard deviation for Gaussian blurring of the saturation channel.
        return_hsv (bool, optional): If True, the intermediate hsv image (masked) is returned.

    Returns:
        list: A list containing prepared images after applying the specified processing steps.
    """
    if isinstance(image, np.ndarray):
        hsv_image = skimage_rgb2hsv(image)
    elif isinstance(image, dask.array.Array):
        hsv_image = dask_rgb2hsv(image)

    if type(image) is dask.array.Array:
        s_image_gaussian = dask_gaussian(hsv_image[..., 1], sigma=mask_gaussian_blur)
    else:
        s_image_gaussian = skimage_gaussian(hsv_image[..., 1], sigma=mask_gaussian_blur)
    
    if type(image) is dask.array.Array:
        thresh = dask_threshold_otsu(s_image_gaussian)
    else:
        thresh = threshold_otsu(s_image_gaussian)
    s_image_gaussian_binary = s_image_gaussian > thresh
    
    s_image_gaussian_binary = binary_fill_holes(s_image_gaussian_binary).astype(int)
    image_out = image * s_image_gaussian_binary[..., np.newaxis]
    hsv_image_out = hsv_image * s_image_gaussian_binary[..., np.newaxis]

    image_out = ((image_out / image_out.max()) * 255).astype(int)
    hsv_image_out = ((hsv_image_out / hsv_image_out.max()) * 255).astype(int)

    if not keep_black_background:
        image_out = np.where(
            image_out == np.array([[0, 0, 0]]),
            np.array([[255, 255, 255]]),
            image_out,
        )
        if return_hsv:
            hsv_image_out = np.where(
                hsv_image_out == np.array([[0, 0, 0]]),
                np.array([[255, 255, 255]]),
                hsv_image_out,
            )

    if return_hsv:
        return image_out, hsv_image_out
    else:
        return image_out
    


def is_grayscale(image):
    """
    Check if an image is grayscale (2D).

    Parameters:
        image (numpy.ndarray): Input image.

    Returns:
        bool: True if the image is grayscale (2D), False otherwise.
    """
    if isinstance(image, np.ndarray):
        if len(image.shape) == 2:
            return True
        elif len(image.shape) == 3 and image.shape[2] == 1:
            return True
    return False