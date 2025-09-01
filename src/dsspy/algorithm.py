import numpy as np
import itertools

from .core import Residue, HBond, HelixPositionType, HelixType, StructureType

# Constants from the C++ implementation
RADIUS_N = 1.65
RADIUS_CA = 1.87
RADIUS_C = 1.76
RADIUS_O = 1.40
RADIUS_SIDE_ATOM = 1.80
RADIUS_WATER = 1.40


def _generate_fibonacci_sphere(n_points: int):
    """
    Generates a sphere of approximately evenly distributed points using the
    Fibonacci sphere algorithm.
    """
    points = np.zeros((n_points, 3))
    phi = np.pi * (3.0 - np.sqrt(5.0))  # Golden angle in radians

    for i in range(n_points):
        y = 1 - (i / float(n_points - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # Golden angle increment

        x = np.cos(theta) * radius
        z = np.sin(theta) * radius

        points[i] = [x, y, z]

    return points


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


MINIMAL_DISTANCE = 0.5
MIN_HBOND_ENERGY = -9.9
COUPLING_CONSTANT = -27.888  # kcal/mol, = -332 * 0.42 * 0.2


def calculate_h_bond_energy(donor: Residue, acceptor: Residue):
    """Calculate the H-bond energy between two residues."""

    energy = 0.0
    if donor.resname != "PRO":
        dist_ho = np.linalg.norm(donor.h_coord - acceptor.o_coord)
        dist_hc = np.linalg.norm(donor.h_coord - acceptor.c_coord)
        dist_nc = np.linalg.norm(donor.n_coord - acceptor.c_coord)
        dist_no = np.linalg.norm(donor.n_coord - acceptor.o_coord)

        if dist_ho < MINIMAL_DISTANCE or dist_hc < MINIMAL_DISTANCE or dist_nc < MINIMAL_DISTANCE or dist_no < MINIMAL_DISTANCE:
            energy = MIN_HBOND_ENERGY
        else:
            energy = COUPLING_CONSTANT / dist_ho - COUPLING_CONSTANT / dist_hc + COUPLING_CONSTANT / dist_nc - COUPLING_CONSTANT / dist_no

        energy = round(energy * 1000) / 1000
        if energy < MIN_HBOND_ENERGY:
            energy = MIN_HBOND_ENERGY

    # This logic matches the C++ source, but the variable names are confusing.
    # In the C++ code, mHBondAcceptor seems to store bonds where the residue is the DONOR.
    if energy < donor.hbond_acceptor[0].energy:
        donor.hbond_acceptor[1] = donor.hbond_acceptor[0]
        donor.hbond_acceptor[0] = HBond(acceptor, energy)
    elif energy < donor.hbond_acceptor[1].energy:
        donor.hbond_acceptor[1] = HBond(acceptor, energy)

    # In the C++ code, mHBondDonor seems to store bonds where the residue is the ACCEPTOR.
    if energy < acceptor.hbond_donor[0].energy:
        acceptor.hbond_donor[1] = acceptor.hbond_donor[0]
        acceptor.hbond_donor[0] = HBond(donor, energy)
    elif energy < acceptor.hbond_donor[1].energy:
        acceptor.hbond_donor[1] = HBond(donor, energy)


def calculate_h_bonds(residues: list[Residue]) -> None:
    """Calculate H-bond energies for all pairs of residues."""

    for res in residues:
        res.assign_hydrogen()

    for i in range(len(residues)):
        for j in range(i + 1, len(residues)):
            if np.linalg.norm(residues[i].ca_coord - residues[j].ca_coord) > 9.0:
                continue

            calculate_h_bond_energy(residues[i], residues[j])
            if j != i + 1:
                calculate_h_bond_energy(residues[j], residues[i])


def _get_atom_spec(residue: Residue):
    """Generator for atom coordinates and their radii."""
    yield residue.n_coord, RADIUS_N
    yield residue.ca_coord, RADIUS_CA
    yield residue.c_coord, RADIUS_C
    yield residue.o_coord, RADIUS_O
    for atom in residue.biopython_residue:
        if atom.get_name() not in ["N", "CA", "C", "O", "H"]:
            # In the C++ code, all side chain atoms that are not H are treated with a generic radius.
            yield atom.get_coord(), RADIUS_SIDE_ATOM

def calculate_accessibility(
    residues: list[Residue],
    n_sphere_points: int = 200,
    water_radius: float = RADIUS_WATER,
):
    """
    Calculates the solvent accessibility for each residue in a list.
    Updates the `accessibility` attribute of each Residue object.
    """
    sphere_points = _generate_fibonacci_sphere(n_sphere_points)

    # Pre-calculate residue centers and radii for neighbor search
    residue_extents = []
    for res in residues:
        all_atom_coords = np.array([spec[0] for spec in _get_atom_spec(res)])
        center = np.mean(all_atom_coords, axis=0)
        radius = np.max(np.linalg.norm(all_atom_coords - center, axis=1))
        residue_extents.append({'center': center, 'radius': radius})

    for i, residue in enumerate(residues):
        # Find neighboring residues
        neighbors = []
        res_i_center = residue_extents[i]['center']
        res_i_radius = residue_extents[i]['radius']

        for j, other_residue in enumerate(residues):
            res_j_center = residue_extents[j]['center']
            res_j_radius = residue_extents[j]['radius']

            # A generous cutoff, similar to the C++ version's logic.
            # Max possible interaction distance is sum of radii + diameter of water.
            cutoff = res_i_radius + res_j_radius + 2 * water_radius

            dist_sq = np.sum((res_i_center - res_j_center)**2)
            if dist_sq < cutoff**2:
                neighbors.append(other_residue)

        # Calculate accessibility for the current residue
        total_accessibility = 0.0
        for atom_coord, atom_radius in _get_atom_spec(residue):
            total_accessibility += _calculate_atom_accessibility(
                atom_coord, atom_radius, neighbors, sphere_points, water_radius
            )
        residue.accessibility = total_accessibility


def _calculate_atom_accessibility(
    atom_coord: np.ndarray,
    atom_radius: float,
    neighboring_residues: list[Residue],
    sphere_points: np.ndarray,
    water_radius: float,
) -> float:
    """Calculates the solvent accessible surface area for a single atom."""

    # Collect all atoms from neighboring residues
    occluding_atoms = list(itertools.chain.from_iterable(_get_atom_spec(res) for res in neighboring_residues))

    # Add water radius to all radii
    radius = atom_radius + water_radius

    # Filter occluding atoms to only those close enough to matter
    # This is a heuristic to speed up the calculation.
    # An atom can't occlude if it's further away than the diameter of the test sphere.
    occluding_atoms_filtered = []
    for occ_coord, occ_radius in occluding_atoms:
        dist_sq = np.sum((atom_coord - occ_coord)**2)
        max_dist = radius + (occ_radius + water_radius)
        if dist_sq < max_dist**2:
             occluding_atoms_filtered.append((occ_coord, occ_radius + water_radius))

    n_points = len(sphere_points)
    accessible_points = 0

    test_points = atom_coord + sphere_points * radius

    for test_point in test_points:
        is_accessible = True
        for occ_coord, occ_radius_with_water in occluding_atoms_filtered:
            dist_sq = np.sum((test_point - occ_coord)**2)
            if dist_sq < occ_radius_with_water**2:
                is_accessible = False
                break
        if is_accessible:
            accessible_points += 1

    # The area of each point is the total sphere area divided by the number of points
    sphere_area = 4 * np.pi * radius**2
    return (accessible_points / n_points) * sphere_area
