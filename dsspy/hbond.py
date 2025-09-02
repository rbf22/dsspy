import numpy as np

from .core import Residue, HBond
from .constants import MINIMAL_DISTANCE, MIN_HBOND_ENERGY, COUPLING_CONSTANT


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
