import numpy as np


def find_fiducial(image: np.ndarray, model_path: str, device: str = "cpu", prob_threshold: float = 0.5) -> np.ndarray:
    """
    Detect fiducial points in an image using a YOLO-based object detection model.

    Args:
        image (np.ndarray): Input image as a numpy array.
        model_path (str): Path to the YOLO model file (e.g., '.pt') for object detection.
        device (str, optional): Device to use for inference ('cpu' or 'cuda'). Default is 'cpu'.
        prob_threshold (float, optional): Probability threshold for filtering detected objects.
          Objects with confidence scores below this threshold will be discarded. Default is 0.5.

    Returns:
        np.ndarray: An array containing the coordinates of detected fiducial points.

    Raises:
        ImportError: If the 'ultralytics' module is not found, an ImportError is raised.

    Note:
        - The 'ultralytics' module is required for this function. You can install it using 'pip install ultralytics'.
        - The function uses the specified YOLO model to detect objects in the input image. Detected objects
          with confidence scores below the 'prob_threshold' are discarded, and the centers of the remaining
          objects are returned as fiducial points.
        - The returned array has shape (N, 2), where N is the number of detected fiducial points.
    """

    try:
        from ultralytics import YOLO
    except ImportError:
        raise ImportError(
            """Could not find module scikit-image.
                          Please run 'pip install ultralytics'"""
        )

    # Create output variables
    passed_matches = []
    centers = []

    model = YOLO(model_path)
    results = model(image, device=device)

    if len(results) == 0:
        return np.array([[0, 0]])

    for i in range(len(results)):
        _r = results[i].boxes.data.cpu().numpy()
        if _r.shape[0] > 0:
            _r_filter = _r[_r[:, 4] > prob_threshold]
            for _r_f in _r_filter:
                if _r_f is not None:
                    passed_matches.append([i] + _r_f.tolist())

    passed_matches = np.array(passed_matches)[:, 1:5]
    for _passed in passed_matches:
        centers.append([((_passed[2] + _passed[0]) / 2), ((_passed[3] + _passed[1]) / 2)])

    return np.array(centers)


def calculate_distance(point1: np.ndarray, point2: np.ndarray) -> float:
    """
    Calculate the Euclidean distance between two 2D points.

    Args:
        point1 (np.ndarray): first list of points
        point2 (np.ndarray): second list of points

    Returns:
        float: distance between the point 1 and 2
    """
    x1, y1 = point1
    x2, y2 = point2
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def correspondences_fiducials(
    points_a: np.ndarray, points_b: np.ndarray, distance_threshold: float = 200
) -> (np.ndarray, np.ndarray):
    """
    Find pairs of points from list1 and list2 that are within the specified threshold distance.

    Args:
        list1: List of 2D coordinates (x, y).
        list2: List of 2D coordinates (x, y).
        threshold: Maximum distance for a pair of points to be considered within the threshold.

    Returns:
        pairs: A list of tuples containing pairs of points that meet the criteria.
    """
    paired_a = []
    paired_b = []

    for point1 in points_a:
        for point2 in points_b:
            distance = calculate_distance(point1, point2)
            if distance <= distance_threshold:
                paired_a.append(point1)
                paired_b.append(point2)

    return np.array(paired_a), np.array(paired_b)
