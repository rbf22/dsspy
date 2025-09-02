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
    Fibonacci sphere algorithm - matching the C++ implementation.
    """
    P = 2 * n_points + 1
    golden_ratio = (1 + np.sqrt(5.0)) / 2
    weight = (4 * np.pi) / P
    
    points = []
    for i in range(-n_points, n_points + 1):
        lat = np.arcsin((2.0 * i) / P)
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


def _atom_intersects_box(atom_coord: np.ndarray, atom_radius: float, box_min: np.ndarray, box_max: np.ndarray) -> bool:
    """Check if an atom intersects with a bounding box."""
    return (atom_coord[0] + atom_radius >= box_min[0] and
            atom_coord[0] - atom_radius <= box_max[0] and
            atom_coord[1] + atom_radius >= box_min[1] and
            atom_coord[1] - atom_radius <= box_max[1] and
            atom_coord[2] + atom_radius >= box_min[2] and
            atom_coord[2] - atom_radius <= box_max[2])


def _calculate_residue_bounding_box(residue: Residue, water_radius: float = RADIUS_WATER):
    """Calculate bounding box for a residue, matching C++ ExtendBox logic."""
    box_min = np.array([np.inf, np.inf, np.inf])
    box_max = np.array([-np.inf, -np.inf, -np.inf])
    
    for atom_coord, atom_radius in _get_atom_spec(residue):
        radius_with_water = atom_radius + 2 * water_radius
        
        box_min = np.minimum(box_min, atom_coord - radius_with_water)
        box_max = np.maximum(box_max, atom_coord + radius_with_water)
    
    return box_min, box_max


class Candidate:
    """Matches the C++ accumulator::candidate structure."""
    def __init__(self, location: np.ndarray, radius_sq: float, distance_sq: float):
        self.location = location
        self.radius_sq = radius_sq  # radius squared for efficiency
        self.distance_sq = distance_sq


def _accumulate_occluding_atoms(atom_coord: np.ndarray, atom_radius: float, 
                               neighboring_residues: list[Residue], 
                               water_radius: float = RADIUS_WATER) -> list[Candidate]:
    """
    Accumulate occluding atoms, matching the C++ accumulator logic.
    """
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
                
                if distance_sq < test_sq and distance_sq > 0.0001:
                    candidate = Candidate(
                        location=occ_coord - atom_coord,
                        radius_sq=r_with_water * r_with_water,
                        distance_sq=distance_sq
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
    """
    Calculates the solvent accessibility for each residue in a list.
    Updates the `accessibility` attribute of each Residue object.
    Now matches C++ implementation more closely.
    """
    sphere_points, weight = _generate_fibonacci_sphere(n_sphere_points)

    for i, residue in enumerate(residues):
        # Calculate accessibility for the current residue
        total_accessibility = 0.0
        for atom_coord, atom_radius in _get_atom_spec(residue):
            total_accessibility += _calculate_atom_accessibility(
                atom_coord, atom_radius, residues, sphere_points, weight, water_radius
            )
        residue.accessibility = total_accessibility


def _calculate_atom_accessibility(
    atom_coord: np.ndarray,
    atom_radius: float,
    all_residues: list[Residue],
    sphere_points: np.ndarray,
    weight: float,
    water_radius: float,
) -> float:
    """Calculates the solvent accessible surface area for a single atom - C++ compatible."""
    
    # Collect occluding atoms using C++ logic
    candidates = _accumulate_occluding_atoms(atom_coord, atom_radius, all_residues, water_radius)
    
    radius = atom_radius + water_radius
    surface = 0.0

    # Generate test points on the sphere surface
    test_points = atom_coord + sphere_points * radius

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
            surface += weight

    # Return surface area (weight already accounts for point density)
    return surface * radius * radius


def calculate_pp_helices(residues: list[Residue], stretch_length: int = 3):
    """
    Identifies Polyproline II (PPII) helices based on phi/psi angles.
    """
    n_residues = len(residues)
    if n_residues < stretch_length:
        return

    epsilon = 29.0
    phi_min = -75.0 - epsilon
    phi_max = -75.0 + epsilon
    psi_min = 145.0 - epsilon
    psi_max = 145.0 + epsilon

    for i in range(n_residues - stretch_length + 1):
        in_helix = True
        for j in range(stretch_length):
            res = residues[i + j]
            if not (res.phi and res.psi and
                    phi_min <= res.phi <= phi_max and
                    psi_min <= res.psi <= psi_max):
                in_helix = False
                break

        if in_helix:
            for j in range(stretch_length):
                res = residues[i + j]
                if res.secondary_structure == StructureType.LOOP:
                    res.secondary_structure = StructureType.HELIX_PPII

                if j == 0: # Start of the helix
                    if res.helix_flags[HelixType.PP] == HelixPositionType.NONE:
                        res.helix_flags[HelixType.PP] = HelixPositionType.START
                    elif res.helix_flags[HelixType.PP] == HelixPositionType.END:
                         res.helix_flags[HelixType.PP] = HelixPositionType.START_AND_END
                elif j == stretch_length - 1: # End of the helix
                    res.helix_flags[HelixType.PP] = HelixPositionType.END
                else: # Middle of the helix
                    res.helix_flags[HelixType.PP] = HelixPositionType.MIDDLE
