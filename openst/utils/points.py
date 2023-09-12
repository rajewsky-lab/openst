from typing import Union

import numpy as np
from scipy.spatial import ConvexHull, Delaunay


def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    if not isinstance(hull, Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p) >= 0


def split_points_by_categorical(points: np.ndarray, cat: np.ndarray) -> dict:
    """
    Split points based on categorical values.

    Args:
        points (np.ndarray): Array of points to be split.
        cat (np.ndarray): Categorical values corresponding to the points.

    Returns:
        dict: A dictionary where keys are unique categorical values and values are arrays of points.

    Raises:
        ValueError: If the provided 'cat' is not 1-dimensional.
        ValueError: If the lengths of 'points' and 'cat' do not match.
    """
    if cat.ndim != 1:
        raise ValueError("The provided 'cat' must be 1-dimensional")
    if len(points) != len(cat):
        raise ValueError("The provided 'points' and 'cat' have different shapes")

    splitted_points = {}

    cat_unique = np.unique(cat)
    for _i_cat_unique in cat_unique:
        splitted_points[_i_cat_unique] = points[cat == _i_cat_unique]

    return splitted_points


def point_inside_which_pointsets(points: np.ndarray, pointsets: Union[list, dict, np.ndarray]):
    """
    Determine which convex hull a set of points lies inside.

    Args:
        points (np.ndarray): Array of points to be checked.
        pointsets (Union[list, dict, np.ndarray]): Convex hulls represented as pointsets.

    Returns:
        np.ndarray: An array indicating which convex hull each point lies inside.

    Raises:
        ValueError: If 'points' is not 2-dimensional.
    """
    if isinstance(pointsets, np.ndarray):
        if pointsets.ndim != 2:
            raise ValueError("'points' must be 2-dimensional")

        pointsets = [pointsets]

    if isinstance(pointsets, list):
        poinsets_keys = np.arange(len(pointsets))
        pointsets = {i: hull for i, hull in enumerate(pointsets)}
    elif isinstance(pointsets, dict):
        poinsets_keys = pointsets.keys()

    hulls = [ConvexHull(pointset) for pointset in pointsets.values()]

    if points.ndim != 2:
        if points.ndim == 1:
            points = points[np.newaxis, :]
        else:
            raise ValueError("Input 'points' must be 2-dimensional")

    points_in_pointsets = np.empty(len(points), dtype="object")

    for _i_p, _p in enumerate(points):
        for poinset_key in poinsets_keys:
            _points_pointsets = pointsets[poinset_key][hulls[poinset_key].vertices]
            if in_hull(_p, _points_pointsets):
                points_in_pointsets[_i_p] = poinset_key
                break

    return points_in_pointsets
