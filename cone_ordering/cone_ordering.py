import math
import heapq
import numpy as np
from typing import List, Tuple
from core.Cone import Cone
from core.distance_matrix import euclidean_distance_between_cones
from core.point_intersection_utils import doIntersect
from config import config

def adjusted_do_intersect(p1, q1, p2, q2, tolerance=1e-3):
    if euclidean_distance_between_cones(p1, q1) < tolerance or euclidean_distance_between_cones(p2, q2) < tolerance:
        return False
    return doIntersect(p1, q1, p2, q2)

def compute_distance_matrix(cones: List[Cone]) -> np.ndarray:
    positions = np.array([(cone.pos.x, cone.pos.y) for cone in cones])
    return np.linalg.norm(positions[:, np.newaxis, :] - positions[np.newaxis, :, :], axis=-1)


def filter_valid_cones(cones: List[Cone], distance_matrix: np.ndarray) -> List[Cone]:
    max_distance = config.MAX_DISTANCE_BETWEEN_CONES
    valid_cones = []
    for i, cone in enumerate(cones):
        distances = distance_matrix[i]
        if any(dist < max_distance for dist in distances):
            valid_cones.append(cone)
    return valid_cones


def compute_neighbouring_cones_matrix(cones: List[Cone], distance_matrix: np.ndarray) -> List[List[Tuple[int, float]]]:
    max_distance = config.MAX_DISTANCE_BETWEEN_CONES  # Ensure this is reasonable
    neighbors = []
    for i, cone in enumerate(cones):
        nearby_cones_and_dists = sorted(
            filter(lambda cd: cd[1] <= max_distance, enumerate(distance_matrix[i])),
            key=lambda cd: cd[1]
        )
        neighbors.append(nearby_cones_and_dists)
    return neighbors


def would_cross_existing_path(tip: Cone, new_cone: Cone, path: List[Cone]) -> bool:
    if len(path) < 2:  # No path to check if fewer than 2 points
        return False
    for p1, p2 in zip(path, path[1:]):
        if adjusted_do_intersect(p1, p2, tip, new_cone, tolerance=8):
            return True
    return False

def compute_curvature(a: Cone, b: Cone, c: Cone) -> float:
    """
    Compute the curvature for three consecutive cones.

    Args:
        a (Cone): The first cone.
        b (Cone): The middle cone.
        c (Cone): The third cone.

    Returns:
        float: The curvature value.
    """
    # Vectors BA and BC
    ba_x, ba_y = a.pos.x - b.pos.x, a.pos.y - b.pos.y
    bc_x, bc_y = c.pos.x - b.pos.x, c.pos.y - b.pos.y

    # Lengths of BA and BC
    ba_len = math.sqrt(ba_x**2 + ba_y**2)
    bc_len = math.sqrt(bc_x**2 + bc_y**2)

    # Handle potential division by zero
    if ba_len == 0 or bc_len == 0:
        return float('inf')

    # Dot product and angle between BA and BC
    dot_product = (ba_x * bc_x + ba_y * bc_y) / (ba_len * bc_len)
    dot_product = max(-1.0, min(1.0, dot_product))  # Clamp to avoid floating-point issues
    angle = math.acos(dot_product)

    # Curvature is the angle deviation from a straight line, normalized by distance
    distance = ba_len + bc_len
    return abs(math.pi - angle) / distance


def compute_angle(a: Cone, b: Cone, c: Cone) -> float:
    """
    Compute the angle (in radians) between the vectors BA and BC.

    Args:
        a (Cone): The first cone.
        b (Cone): The middle cone.
        c (Cone): The third cone.

    Returns:
        float: The angle in radians.
    """
    # Vectors BA and BC
    ba_x, ba_y = a.pos.x - b.pos.x, a.pos.y - b.pos.y
    bc_x, bc_y = c.pos.x - b.pos.x, c.pos.y - b.pos.y

    # Lengths of BA and BC
    ba_len = math.sqrt(ba_x**2 + ba_y**2)
    bc_len = math.sqrt(bc_x**2 + bc_y**2)

    # Handle potential division by zero
    if ba_len == 0 or bc_len == 0:
        return 0  # Treat as invalid (angle of 0)

    # Dot product and angle between BA and BC
    dot_product = (ba_x * bc_x + ba_y * bc_y) / (ba_len * bc_len)
    dot_product = max(-1.0, min(1.0, dot_product))  # Clamp to avoid floating-point issues
    angle = math.acos(dot_product)
    return angle


def get_ordered_list_of_cones(cones: List[Cone], start_cone: Cone) -> List[Cone]:
    """
    Orders a list of cones starting from the specified `start_cone`.

    Args:
        cones (List[Cone]): The list of cones to order.
        start_cone (Cone): The starting cone for this set of cones.

    Returns:
        List[Cone]: The ordered list of cones.
    """
    # Compute the distance matrix for all cones
    distance_matrix = compute_distance_matrix(cones)

    # Filter invalid cones based on distance criteria
    valid_cones = filter_valid_cones(cones, distance_matrix)

    visited = {start_cone}
    path = [start_cone]

    # Priority queue to explore neighbors
    neighbors = compute_neighbouring_cones_matrix(valid_cones, distance_matrix)
    priority_queue = []
    for neighbor_idx, dist in neighbors[valid_cones.index(start_cone)]:
        heapq.heappush(priority_queue, (dist, neighbor_idx, valid_cones[neighbor_idx]))

    # Process the priority queue
    while priority_queue:
        _, _, current_cone = heapq.heappop(priority_queue)
        if current_cone in visited:
            continue

        # Validate angle if path has at least two previous cones
        if len(path) >= 2:
            a, b = path[-2], path[-1]
            angle = compute_angle(a, b, current_cone)
            min_angle = math.radians(30)  # Minimum allowed angle (e.g., 30 degrees)
            if angle < min_angle:
                continue  
            
        # Validate and add cone
        if not would_cross_existing_path(path[-1], current_cone, path):
            path.append(current_cone)
            visited.add(current_cone)

            # Add current cone's neighbors to the queue
            for neighbor_idx, dist in neighbors[valid_cones.index(current_cone)]:
                neighbor_cone = valid_cones[neighbor_idx]
                if neighbor_cone not in visited:
                    heapq.heappush(priority_queue, (dist, neighbor_idx, neighbor_cone))

    # Ensure path loops back to the start if possible
    if len(path) > 1 and euclidean_distance_between_cones(path[-1], path[0]) < config.MAX_DISTANCE_BETWEEN_CONES:
        path.append(path[0])

    return path


def order_blue_and_yellow_cones(blue_cones: List[Cone], yellow_cones: List[Cone]) -> Tuple[List[Cone], List[Cone]]:
    """
    Orders blue and yellow cones into separate tracks, considering angle constraints.

    Args:
        blue_cones (List[Cone]): The list of blue cones.
        yellow_cones (List[Cone]): The list of yellow cones.

    Returns:
        Tuple[List[Cone], List[Cone]]: Ordered lists for blue and yellow cones.
    """
    # Find the closest cone to the origin for each color
    start_blue = min(blue_cones, key=lambda cone: euclidean_distance_between_cones(cone, config.ORIGIN))
    start_yellow = min(yellow_cones, key=lambda cone: euclidean_distance_between_cones(cone, config.ORIGIN))

    # Debug: Print starting cones
    print(f"Start Blue: ({start_blue.pos.x}, {start_blue.pos.y})")
    print(f"Start Yellow: ({start_yellow.pos.x}, {start_yellow.pos.y})")

    # Order cones separately for blue and yellow
    ordered_blue = get_ordered_list_of_cones(blue_cones, start_blue)
    ordered_yellow = get_ordered_list_of_cones(yellow_cones, start_yellow)

    return ordered_blue, ordered_yellow