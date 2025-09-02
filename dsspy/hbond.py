import numpy as np

from .core import Residue, HBond
from .constants import MINIMAL_DISTANCE, MIN_HBOND_ENERGY, COUPLING_CONSTANT


def calculate_h_bond_energy(donor: Residue, acceptor: Residue):
    """Calculate the H-bond energy between two residues."""
    
    energy = 0.0
    
    # Critical fix: Only calculate energy if donor is NOT proline
    # Proline cannot donate hydrogen bonds because its nitrogen is part of a ring
    if donor.resname != "PRO":
        dist_ho = np.linalg.norm(donor.h_coord - acceptor.o_coord)
        dist_hc = np.linalg.norm(donor.h_coord - acceptor.c_coord)
        dist_nc = np.linalg.norm(donor.n_coord - acceptor.c_coord)
        dist_no = np.linalg.norm(donor.n_coord - acceptor.o_coord)

        if (dist_ho < MINIMAL_DISTANCE or 
            dist_hc < MINIMAL_DISTANCE or 
            dist_nc < MINIMAL_DISTANCE or 
            dist_no < MINIMAL_DISTANCE):
            energy = MIN_HBOND_ENERGY
        else:
            energy = (COUPLING_CONSTANT / dist_ho - 
                     COUPLING_CONSTANT / dist_hc + 
                     COUPLING_CONSTANT / dist_nc - 
                     COUPLING_CONSTANT / dist_no)

        # Round to match DSSP compatibility mode
        energy = round(energy * 1000) / 1000
        
        if energy < MIN_HBOND_ENERGY:
            energy = MIN_HBOND_ENERGY

    # Only update hydrogen bond arrays if energy is significant
    # This prevents weak/zero bonds from overwriting stronger ones
    
    # Update donor's acceptor bonds (bonds where this residue donates)
    if energy < donor.hbond_acceptor[0].energy:
        donor.hbond_acceptor[1] = donor.hbond_acceptor[0]
        donor.hbond_acceptor[0] = HBond(acceptor, energy)
    elif energy < donor.hbond_acceptor[1].energy:
        donor.hbond_acceptor[1] = HBond(acceptor, energy)

    # Update acceptor's donor bonds (bonds where this residue accepts)
    if energy < acceptor.hbond_donor[0].energy:
        acceptor.hbond_donor[1] = acceptor.hbond_donor[0]
        acceptor.hbond_donor[0] = HBond(donor, energy)
    elif energy < acceptor.hbond_donor[1].energy:
        acceptor.hbond_donor[1] = HBond(donor, energy)

    return energy

def assign_hydrogen_to_residues(residues: list[Residue]):
    """
    Assign hydrogen positions for all residues.
    This is critical for accurate H-bond calculation.
    """
    for i, residue in enumerate(residues):
        # Start with nitrogen position
        residue.h_coord = residue.n_coord.copy()
        
        # For non-proline residues with a previous residue
        if residue.resname != "PRO" and i > 0:
            prev_residue = residues[i - 1]
            prev_c = prev_residue.c_coord
            prev_o = prev_residue.o_coord
            
            # Calculate CO vector and normalize it
            co_vector = prev_c - prev_o
            co_distance = np.linalg.norm(co_vector)
            
            if co_distance > 0:
                co_unit = co_vector / co_distance
                # Place hydrogen along the CO vector direction from nitrogen
                residue.h_coord += co_unit

def calculate_h_bonds(residues: list[Residue]) -> None:
    """Calculate H-bond energies for all pairs of residues."""

    # Note: H-bond arrays are already initialized in Residue.__init__()
    # But we need to reset them to ensure clean state
    for res in residues:
        res.hbond_acceptor = [HBond(None, 0.0), HBond(None, 0.0)]
        res.hbond_donor = [HBond(None, 0.0), HBond(None, 0.0)]

    for res in residues:
        res.assign_hydrogen()

    for i in range(len(residues)):
        for j in range(i + 1, len(residues)):
            if np.linalg.norm(residues[i].ca_coord - residues[j].ca_coord) > 9.0:
                continue

            calculate_h_bond_energy(residues[i], residues[j])
            if j != i + 1:
                calculate_h_bond_energy(residues[j], residues[i])
