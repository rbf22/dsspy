from dsspy.core import StructureType, HelixType, HelixPositionType
from dsspy.algorithm import calculate_pp_helices

class MockResidue:
    def __init__(self, phi, psi):
        self.phi = phi
        self.psi = psi
        self.secondary_structure = StructureType.LOOP
        self.helix_flags = {helix_type: HelixPositionType.NONE for helix_type in HelixType}

def test_calculate_pp_helices():
    """
    Tests the calculate_pp_helices function with the C++ ported logic.
    """
    # Create a list of mock residues with phi/psi angles for a PPII helix
    residues = [MockResidue(-75, 145) for _ in range(7)]
    residues[0].phi = 0 # Not in helix
    residues[6].phi = 0 # Not in helix

    calculate_pp_helices(residues, stretch_length=3)

    # Expected state based on the quirky C++ loop (for i in 1..n-4):
    # N, S, M, M, M, E, N
    assert residues[0].helix_flags[HelixType.PP] == HelixPositionType.NONE
    assert residues[1].helix_flags[HelixType.PP] == HelixPositionType.START
    assert residues[2].helix_flags[HelixType.PP] == HelixPositionType.MIDDLE
    assert residues[3].helix_flags[HelixType.PP] == HelixPositionType.MIDDLE
    assert residues[4].helix_flags[HelixType.PP] == HelixPositionType.MIDDLE
    assert residues[5].helix_flags[HelixType.PP] == HelixPositionType.END
    assert residues[6].helix_flags[HelixType.PP] == HelixPositionType.NONE

    assert residues[0].secondary_structure == StructureType.LOOP
    assert residues[1].secondary_structure == StructureType.HELIX_PPII
    assert residues[2].secondary_structure == StructureType.HELIX_PPII
    assert residues[3].secondary_structure == StructureType.HELIX_PPII
    assert residues[4].secondary_structure == StructureType.HELIX_PPII
    assert residues[5].secondary_structure == StructureType.HELIX_PPII
    assert residues[6].secondary_structure == StructureType.LOOP
