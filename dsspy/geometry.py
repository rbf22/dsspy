"""
Geometric calculations for DSSP.
"""

import numpy as np


def _generate_fibonacci_sphere(n_points: int):
    """
    Generates a sphere of approximately evenly distributed points using the
    Fibonacci sphere algorithm - matching the C++ implementation.
    """
    p_value = 2 * n_points + 1
    golden_ratio = (1 + np.sqrt(5.0)) / 2
    weight = (4 * np.pi) / p_value

    points = []
    for i in range(-n_points, n_points + 1):
        lat = np.arcsin((2.0 * i) / p_value)
        lon = (i % golden_ratio) * 2 * np.pi / golden_ratio

        points.append([
            np.sin(lon) * np.cos(lat),
            np.cos(lon) * np.cos(lat),
            np.sin(lat)
        ])

    return np.array(points), weight


def dihedral_angle(p1, p2, p3, p4):
    """
    Calculates the dihedral angle between four points.
    The angle is calculated in degrees.
    """
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    # Normalize b2
    b2_norm = b2 / np.linalg.norm(b2)

    y = np.dot(np.cross(n1, n2), b2_norm)
    x = np.dot(n1, n2)

    angle_rad = np.arctan2(y, x)
    angle_deg = np.rad2deg(angle_rad)
    if np.isclose(angle_deg, 180.0):
        return -180.0
    return angle_deg
