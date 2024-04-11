import numpy as np
from skimage.transform import SimilarityTransform


def apply_transform(in_coords: np.ndarray, transform: SimilarityTransform, check_bounds=False):
    # Check if transform within the acceptable bounds
    if (
        check_bounds
        and ((transform.rotation > np.pi / 4)
        or (transform.scale > 2 or transform.scale < 0.5)
        or (transform.translation.max() > in_coords.max(axis=0).max()))
    ):
        return in_coords

    # If the previous filter passes, apply the transformation
    out_coords = np.dot(
        transform.params,
        np.concatenate(
            [
                in_coords[:, ::-1],
                np.ones((len(in_coords), 1)),
            ],
            axis=1,
        ).T,
    ).T

    return out_coords
