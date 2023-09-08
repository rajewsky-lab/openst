import numpy as np

from skimage.feature import SIFT, match_descriptors

def _find_matches_loftr(im_0: np.ndarray, im_1: np.ndarray, pretrained: str = "outdoor") -> tuple:
    """
    Find matching keypoints between two images using LoFTR.

    Args:
        im_0 (np.ndarray): First input image.
        im_1 (np.ndarray): Second input image.
        pretrained (str, optional): Pretrained LoFTR model variant. Default is 'outdoor'.

    Returns:
        tuple: A tuple containing two arrays:
            - keypoints0 (np.ndarray): Array of keypoints from the first image.
            - keypoints1 (np.ndarray): Array of keypoints from the second image.

    Notes:
        - This function uses LoFTR to find matching keypoints between two input images.
        - The 'pretrained' parameter specifies the variant of the LoFTR model to use.
    """

    try:
        import kornia.feature as KF
        import torch
    except ImportError:
        raise ImportError(
            """Optional modules need to be installed to run this feature.
               Please run 'pip install kornia'"""
        )

    matcher = KF.LoFTR(pretrained=pretrained)

    input_dict = {
        "image0": torch.tensor(im_0)[None, None].float(),
        "image1": torch.tensor(im_1)[None, None].float(),
    }

    with torch.inference_mode():
        correspondences = matcher(input_dict)

    return correspondences["keypoints0"].cpu().numpy(), correspondences["keypoints1"].cpu().numpy()


def _find_matches_sift(im_0: np.ndarray, im_1: np.ndarray) -> tuple:
    """
    Find matching keypoints between two images using SIFT.

    Args:
        im_0 (np.ndarray): First input image.
        im_1 (np.ndarray): Second input image.

    Returns:
        tuple: A tuple containing two arrays:
            - keypoints0 (np.ndarray): Array of keypoints from the first image.
            - keypoints1 (np.ndarray): Array of keypoints from the second image.

    Notes:
        - This function uses SIFT to find matching keypoints between two input images.
    """
    feat_descriptor = SIFT()

    feat_descriptor.detect_and_extract(im_0)
    keypoints0 = feat_descriptor.keypoints
    descriptors0 = feat_descriptor.descriptors

    feat_descriptor.detect_and_extract(im_1)
    keypoints1 = feat_descriptor.keypoints
    descriptors1 = feat_descriptor.descriptors

    matches01 = match_descriptors(descriptors0, descriptors1, max_ratio=0.6, cross_check=True)

    return keypoints0[matches01], keypoints1[matches01]


def find_matches(src, dst, method="LoFTR") -> tuple:
    """
    Find matching keypoints between source and destination images using a specified method.

    Args:
        src (list): List of source images for keypoint matching.
        dst (list): List of destination images for keypoint matching.
        method (str, optional): Method (LoFTR or SIFT) used for feature detection and matching.

    Returns:
        tuple: A tuple containing two arrays:
            - mkpts0 (np.ndarray): Array of keypoints from the destination images.
            - mkpts1 (np.ndarray): Array of keypoints from the source images.

    Raises:
        TypeError: If 'src' is not a list of images or 'dst' is not a list of images.

    Notes:
        - This function uses LoFTR to find matching keypoints between source and destination images.
        - 'src' and 'dst' should be lists of images to be matched.
        - The function returns a tuple of arrays containing the keypoints for matching.
    """

    if type(src) is not list:
        raise TypeError("'src' must be a list of images")

    if type(dst) is not list:
        raise TypeError("'dst' must be a list of images")

    mkpts0 = np.array([[0, 0]])
    mkpts1 = np.array([[0, 0]])

    for _i_dst in dst:
        for _i_src in src:
            if method == "LoFTR":
                _i_mkpts0, _i_mkpts1 = _find_matches_loftr(_i_dst, _i_src)
            elif method == "SIFT":
                _i_mkpts0, _i_mkpts1 = _find_matches_sift(_i_dst, _i_src)
            else:
                raise ValueError(f"Registration method {method} not supported")

            mkpts0 = np.concatenate([_i_mkpts0, mkpts0])
            mkpts1 = np.concatenate([_i_mkpts1, mkpts1])

    return mkpts0, mkpts1