import numpy as np
from enum import Enum

class StructureType(Enum):
    LOOP = ' '
    ALPHA_HELIX = 'H'
    BETA_BRIDGE = 'B'
    STRAND = 'E'
    HELIX_3 = 'G'
    HELIX_5 = 'I'
    HELIX_PPII = 'P'
    TURN = 'T'
    BEND = 'S'

class HelixType(Enum):
    _3_10 = 0
    ALPHA = 1
    PI = 2
    PP = 3

class HelixPositionType(Enum):
    NONE = 0
    START = 1
    MIDDLE = 2
    END = 3
    START_AND_END = 4

class ChainBreakType(Enum):
    NONE = 0
    NEW_CHAIN = 1
    GAP = 2


class HBond:
    def __init__(self, residue, energy):
        self.residue = residue
        self.energy = energy

    def __repr__(self):
        if self.residue:
            return f"<HBond to {self.residue.id} with energy {self.energy}>"
        return f"<HBond to None with energy {self.energy}>"


class Residue:
    def __init__(self, biopython_residue, number):
        self.biopython_residue = biopython_residue
        self.number = number
        self.secondary_structure = StructureType.LOOP
        self.hbond_acceptor = [HBond(None, 0), HBond(None, 0)]
        self.hbond_donor = [HBond(None, 0), HBond(None, 0)]
        self.beta_partner = [None, None]
        self.sheet = 0
        self.strand = 0
        self.helix_flags = {helix_type: HelixPositionType.NONE for helix_type in HelixType}
        self.bend = False
        self.chain_break = ChainBreakType.NONE
        self.alpha = 360.0
        self.kappa = 360.0
        self.phi = 360.0
        self.psi = 360.0
        self.tco = 0.0
        self.omega = 360.0
        self.accessibility = 0.0
        self.prev_residue = None
        self.next_residue = None
        self.h_coord = None

    def assign_hydrogen(self):
        """
        Assigns the coordinates of the hydrogen atom based on the previous residue.
        """
        if self.resname != 'PRO' and self.prev_residue:
            n = self.n_coord
            pc = self.prev_residue.c_coord
            po = self.prev_residue.o_coord

            co_dist = pc - po
            self.h_coord = n + co_dist / np.linalg.norm(co_dist)
        else:
            # For prolines and the first residue in a chain, we don't have a previous C=O to base the H on.
            # The C++ code initializes H to N's position in this case.
            self.h_coord = self.n_coord

    @property
    def ca_coord(self):
        return self.biopython_residue['CA'].get_coord()

    @property
    def n_coord(self):
        return self.biopython_residue['N'].get_coord()

    @property
    def c_coord(self):
        return self.biopython_residue['C'].get_coord()

    @property
    def o_coord(self):
        return self.biopython_residue['O'].get_coord()

    @property
    def id(self):
        return self.biopython_residue.get_id()

    @property
    def resname(self):
        return self.biopython_residue.get_resname()

    def __repr__(self):
        return f"<Residue {self.resname} id={self.id}>"

class Bridge:
    def __init__(self, r1, r2, type):
        self.r1 = r1
        self.r2 = r2
        self.type = type
        self.ladder = None
        self.sheet = None

    def __repr__(self):
        return f"<Bridge {self.r1.id}-{self.r2.id} type={self.type}>"
