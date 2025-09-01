import numpy as np

from dsspy.core import Residue, HBond


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
