"""
This module contains functions for calculating solvent accessibility of residues.
"""
from dataclasses import dataclass
import numpy as np

from .core import Residue
from .constants import (
    RADIUS_N,
    RADIUS_CA,
    RADIUS_C,
    RADIUS_O,
    RADIUS_SIDE_ATOM,
    RADIUS_WATER,
)
from .geometry import _generate_fibonacci_sphere


def _get_atom_spec(residue: Residue):
    """Generator for atom coordinates and their radii."""
    yield residue.n_coord, RADIUS_N
    yield residue.ca_coord, RADIUS_CA
    yield residue.c_coord, RADIUS_C
    yield residue.o_coord, RADIUS_O
    for atom in residue.biopython_residue:
        if atom.get_name() not in ["N", "CA", "C", "O", "H"]:
            # In the C++ code, all side chain atoms that are not H are treated
            # with a generic radius.
            yield atom.get_coord(), RADIUS_SIDE_ATOM


def _atom_intersects_box(
    atom_coord: np.ndarray, atom_radius: float, box_min: np.ndarray, box_max: np.ndarray
) -> bool:
    """Check if an atom intersects with a bounding box."""
    return (
        box_min[0] <= atom_coord[0] + atom_radius
        and box_max[0] >= atom_coord[0] - atom_radius
        and box_min[1] <= atom_coord[1] + atom_radius
        and box_max[1] >= atom_coord[1] - atom_radius
        and box_min[2] <= atom_coord[2] + atom_radius
        and box_max[2] >= atom_coord[2] - atom_radius
    )


def _calculate_residue_bounding_box(
    residue: Residue, water_radius: float = RADIUS_WATER
):
    """Calculate bounding box for a residue, matching C++ ExtendBox logic."""
    box_min = np.array([np.inf, np.inf, np.inf])
    box_max = np.array([-np.inf, -np.inf, -np.inf])

    for atom_coord, atom_radius in _get_atom_spec(residue):
        radius_with_water = atom_radius + 2 * water_radius

        box_min = np.minimum(box_min, atom_coord - radius_with_water)
        box_max = np.maximum(box_max, atom_coord + radius_with_water)

    return box_min, box_max


@dataclass
class Sphere:
    """Represents a sphere of points for accessibility calculation."""
    points: np.ndarray
    weight: float


class Candidate:
    """Matches the C++ accumulator::candidate structure."""
    # pylint: disable=too-few-public-methods
    def __init__(self, location: np.ndarray, radius_sq: float, distance_sq: float):
        self.location = location
        self.radius_sq = radius_sq  # radius squared for efficiency
        self.distance_sq = distance_sq


def _accumulate_occluding_atoms(
    atom_coord: np.ndarray,
    atom_radius: float,
    neighboring_residues: list[Residue],
    water_radius: float = RADIUS_WATER,
) -> list[Candidate]:
    """
    Accumulate occluding atoms, matching the C++ accumulator logic.
    """
    # pylint: disable=too-many-locals
    candidates = []
    d_with_water = atom_radius + water_radius

    for residue in neighboring_residues:
        # Calculate residue bounding box
        box_min, box_max = _calculate_residue_bounding_box(residue, water_radius)

        # Check if atom intersects with residue bounding box
        if _atom_intersects_box(atom_coord, d_with_water, box_min, box_max):
            for occ_coord, occ_radius in _get_atom_spec(residue):
                r_with_water = occ_radius + water_radius

                distance_sq = np.sum((atom_coord - occ_coord) ** 2)

                test_radius = d_with_water + r_with_water
                test_sq = test_radius * test_radius

                if 0.0001 < distance_sq < test_sq:
                    candidate = Candidate(
                        location=occ_coord - atom_coord,
                        radius_sq=r_with_water * r_with_water,
                        distance_sq=distance_sq,
                    )
                    candidates.append(candidate)

    # Sort by distance (matching C++ sort_heap behavior)
    candidates.sort(key=lambda x: x.distance_sq)
    return candidates


def calculate_accessibility(
    residues: list[Residue],
    n_sphere_points: int = 200,
    water_radius: float = RADIUS_WATER,
):
    """Calculates and assigns solvent accessibility for each residue.

    This function implements a Shrake & Rupley style algorithm to calculate
    the solvent accessible surface area (SASA) for each residue. It does this
    by creating a sphere of points around each atom and checking how many of
    these points are occluded by other atoms.

    Note: This function has side effects, as it modifies the `accessibility`
    attribute on the Residue objects in the input list.

    Args:
        residues (list[Residue]): A list of all Residue objects in the structure.
        n_sphere_points (int, optional): The number of points to generate on the
            surface of the sphere around each atom for the accessibility check.
            More points are more accurate but slower. Defaults to 200.
        water_radius (float, optional): The radius of a water molecule, used as
            the probe radius. Defaults to RADIUS_WATER.
    """
    sphere_points, weight = _generate_fibonacci_sphere(n_sphere_points)
    sphere = Sphere(points=sphere_points, weight=weight)

    for residue in residues:
        # Calculate accessibility for the current residue
        total_accessibility = 0.0
        for atom_coord, atom_radius in _get_atom_spec(residue):
            total_accessibility += _calculate_atom_accessibility(
                atom_coord,
                atom_radius,
                residues,
                sphere,
                water_radius,
            )
        residue.accessibility = total_accessibility


def _calculate_atom_accessibility(
    atom_coord: np.ndarray,
    atom_radius: float,
    all_residues: list[Residue],
    sphere: Sphere,
    water_radius: float,
) -> float:
    """Calculates the solvent accessible surface area for a single atom - C++ compatible."""

    # Collect occluding atoms using C++ logic
    candidates = _accumulate_occluding_atoms(atom_coord, atom_radius, all_residues, water_radius)

    radius = atom_radius + water_radius
    surface = 0.0

    # Generate test points on the sphere surface
    test_points = atom_coord + sphere.points * radius

    for test_point in test_points:
        is_accessible = True

        # Check against all occluding candidates
        for candidate in candidates:
            # Calculate distance squared from test point to occluding atom center
            dist_sq = np.sum((test_point - (atom_coord + candidate.location)) ** 2)

            if dist_sq < candidate.radius_sq:
                is_accessible = False
                break

        if is_accessible:
            surface += sphere.weight

    # Return surface area (weight already accounts for point density)
    return surface * radius * radius
