import logging
from collections.abc import Callable
from itertools import product

import numpy as np
from skimage.measure import ransac
from skimage.transform import SimilarityTransform

SUPPORTED_MATCHING_METHODS = ["LoFTR", "SIFT", "KeyNet"]


def _find_matches_loftr(im_0: np.ndarray, im_1: np.ndarray, pretrained: str = "outdoor", device: str = "cpu") -> tuple:
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
            """Optional modules need to be installed to run thi s feature.
               Please run 'pip install kornia'"""
        )

    matcher = KF.LoFTR(pretrained=pretrained).to(device)

    input_dict = {
        "image0": torch.tensor(im_0.copy())[None, None].float().to(device),
        "image1": torch.tensor(im_1.copy())[None, None].float().to(device),
    }

    with torch.inference_mode():
        correspondences = matcher(input_dict)

    return correspondences["keypoints0"].cpu().numpy(), correspondences["keypoints1"].cpu().numpy()


def _find_matches_keynet(im_0: np.ndarray, im_1: np.ndarray, device: str = "cpu") -> tuple:
    """
    Find matching keypoints between two images using LoFTR.

    Args:
        im_0 (np.ndarray): First input image.
        im_1 (np.ndarray): Second input image.

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
            """Optional modules need to be installed to run thi s feature.
               Please run 'pip install kornia'"""
        )

    feature = KF.KeyNetAffNetHardNet(5000, True).eval().to(device)

    input_dict = {
        "image0": torch.tensor(im_0)[None, None].float().to(device),
        "image1": torch.tensor(im_1)[None, None].float().to(device),
    }

    hw1 = torch.tensor(input_dict['image0'].shape[2:])
    hw2 = torch.tensor(input_dict['image0'].shape[2:])

    adalam_config = {"device": device}

    with torch.inference_mode():
        lafs1, resps1, descs1 = feature(input_dict['image0'])
        lafs2, resps2, descs2 = feature(input_dict['image1'])
        dists, idxs = KF.match_adalam(
            descs1.squeeze(0),
            descs2.squeeze(0),
            lafs1,
            lafs2,  # Adalam takes into account also geometric information
            config=adalam_config,
            hw1=hw1,
            hw2=hw2,  # Adalam also benefits from knowing image size
        )

    def get_matching_keypoints(lafs1, lafs2, idxs):
        mkpts1 = KF.get_laf_center(lafs1).squeeze()[idxs[:, 0]].detach().cpu().numpy()
        mkpts2 = KF.get_laf_center(lafs2).squeeze()[idxs[:, 1]].detach().cpu().numpy()

        return mkpts1, mkpts2

    mkpts1, mkpts2 = get_matching_keypoints(lafs1, lafs2, idxs)

    return mkpts1, mkpts2


def _find_matches_sift(im_0: np.ndarray, im_1: np.ndarray, *args, **kwargs) -> tuple:
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
    try:
        from skimage.feature import SIFT, match_descriptors
    except ImportError:
        raise ImportError(
            """Could not find module scikit-image.
                          Please run 'pip install scikit-image'"""
        )
    feat_descriptor = SIFT()

    feat_descriptor.detect_and_extract(im_0)
    keypoints0 = feat_descriptor.keypoints
    descriptors0 = feat_descriptor.descriptors

    feat_descriptor.detect_and_extract(im_1)
    keypoints1 = feat_descriptor.keypoints
    descriptors1 = feat_descriptor.descriptors

    matches01 = match_descriptors(descriptors0, descriptors1, max_ratio=0.6, cross_check=True)

    return keypoints0[matches01[:, 0]], keypoints1[matches01[:, 1]]


def find_matches(
    src: list,
    dst: list,
    method="LoFTR",
    prefilter: bool = False,
    ransac_min_samples: float = 1,
    ransac_residual_threshold: float = 1,
    ransac_max_trials: float = 10000,
    device: str = "cpu",
) -> tuple:
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

    if method not in SUPPORTED_MATCHING_METHODS:
        raise ValueError(f"Feature matching method '{method}' is not supported")

    mkpts0 = np.array([])
    mkpts1 = np.array([])

    for _i_dst in dst:
        for _i_src in src:
            if method == "LoFTR":
                _i_mkpts0, _i_mkpts1 = _find_matches_loftr(_i_dst, _i_src, device=device)
            elif method == "SIFT":
                _i_mkpts0, _i_mkpts1 = _find_matches_sift(_i_dst, _i_src)
            elif method == 'KeyNet':
                _i_mkpts0, _i_mkpts1 = _find_matches_keynet(_i_dst, _i_src, device=device)
            else:
                raise ValueError(f"Registration method {method} not supported")

            if prefilter and _i_mkpts0.size > 0:
                _, inliers = ransac(
                    (_i_mkpts0, _i_mkpts1),
                    SimilarityTransform,
                    min_samples=ransac_min_samples,
                    residual_threshold=ransac_residual_threshold,
                    max_trials=ransac_max_trials,
                )

                if inliers is not None:
                    inliers = inliers > 0
                    _i_mkpts0 = _i_mkpts0[inliers.flatten()]
                    _i_mkpts1 = _i_mkpts1[inliers.flatten()]
                else:
                    _i_mkpts0 = []
                    _i_mkpts1 = []

            if len(_i_mkpts0) > 0 and len(mkpts0) > 0:
                mkpts0 = np.concatenate([_i_mkpts0, mkpts0])
                mkpts1 = np.concatenate([_i_mkpts1, mkpts1])
            elif len(_i_mkpts0) > 0 and len(mkpts0) == 0:
                mkpts0 = _i_mkpts0
                mkpts1 = _i_mkpts1

    return mkpts0, mkpts1


def match_images(
    src: np.ndarray,
    dst: np.ndarray,
    feature_matcher: str = "LoFTR",
    flips: list = [[1, 1], [1, -1], [-1, 1], [-1, -1]],
    rotations: list = [0, 90],
    src_augmenter: Callable = None,
    dst_augmenter: Callable = None,
    prefilter: bool = False,
    ransac_enabled: bool = True,
    ransac_min_samples: float = 1,
    ransac_residual_threshold: float = 1,
    ransac_max_trials: float = 10000,
    device: str = "cpu",
) -> (np.ndarray, np.ndarray, list, float):
    """
    Matching of two images (A,B), with augmentation (optionally)

    Args:
        src (np.ndarray): Source image.
        dst (np.ndarray): Destination image.
        feature_matcher (str): Method used for feature detection and matching.
        flips (list): a 2d list of flips (vertical, horizontal) applied to dst before matching.
        rotations (list): a list of rotations applied to dst before matching.
        src_augmenter (function): augmentation function to apply to image A.
        dst_augmenter (function): augmentation function to apply to image B.

    Returns:
        tuple: A tuple containing:
            - The (matching) keypoints of image A, at the best combination of flips/rotations.
            - The (matching) keypoints of image B.
            - The best flip of image B (yielding highest A/B matches)
            - The best rotation of image B (yielding highest A/B matches)
    """
    # TODO: check flips list
    # TODO: check rotations list
    # TODO: documentation

    if (isinstance(src, list) and isinstance(src_augmenter, Callable)) \
        or (isinstance(dst, list) and isinstance(dst_augmenter, Callable)):
        raise NotImplementedError("Cannot run matching >1 src/dst images and an augmenter\n" +
                                  "Please choose only one augmentation strategy")

    _src, _dst = src, dst
    max_keypoints = 0
    best_flip = flips[0]
    best_rotation = rotations[0]
    _best_mkpts0 = None
    _best_mkpts1 = None

    for (_flip_x, _flip_y), _rotation in product(flips, rotations):
        logging.info(f"Flip {_flip_x, _flip_y}, rotation {_rotation}")

        # Preparing image and pseudoimage modalities for the feature matching model
        if callable(src_augmenter):
            _src = src_augmenter(src, flip=[_flip_x, _flip_y], rotation=_rotation)
        if callable(dst_augmenter):
            _dst = dst_augmenter(dst, flip=[_flip_x, _flip_y], rotation=_rotation)

        # Find matching keypoints between image and STS pseudoimage modalities
        mkpts0, mkpts1 = find_matches(_src, _dst, feature_matcher, prefilter, device=device)

        logging.info(f"{len(mkpts0)} matches")

        # Run RANSAC to remove outliers,
        # after all pairwise channels have been matched
        if ransac_enabled:
            sum_inliers = 0
            best_inliers = None
            for _ in range(ransac_max_trials):
                _, inliers = ransac(
                    (mkpts0, mkpts1),
                    SimilarityTransform,
                    min_samples=ransac_min_samples,
                    residual_threshold=ransac_residual_threshold,
                    max_trials=1000,
                )
                inliers = inliers > 0

                if inliers.sum() > sum_inliers:
                    sum_inliers = inliers.sum()
                    best_inliers = inliers
        else:
            best_inliers = np.ones(len(mkpts0), dtype=bool)

        if best_inliers.sum() > max_keypoints:
            max_keypoints = best_inliers.sum()
            best_flip = [_flip_x, _flip_y]
            best_rotation = _rotation

            _best_mkpts0 = mkpts0[best_inliers.flatten()]
            _best_mkpts1 = mkpts1[best_inliers.flatten()]

        logging.info(f"{best_inliers.sum()} inliers (RANSAC)")

    # Retrieve the results for the best flip combination
    # Filter keypoints with selected inliers
    in_mkpts0 = _best_mkpts0
    in_mkpts1 = _best_mkpts1

    return in_mkpts0, in_mkpts1, best_flip, best_rotation
