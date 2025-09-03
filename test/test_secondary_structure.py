"""
Tests for the secondary_structure module.
"""
from unittest.mock import Mock
from dsspy.core import StructureType, HelixType, HelixPositionType
from dsspy.secondary_structure import calculate_pp_helices

def create_mock_residue(phi=None, psi=None):
    """Creates a mock residue for testing."""
    res = Mock()
    res.phi = phi
    res.psi = psi
    res.secondary_structure = StructureType.LOOP
    res.helix_flags = {ht: HelixPositionType.NONE for ht in HelixType}
    return res

def test_calculate_pp_helices():
    """
    Tests the calculate_pp_helices function.
    """
    # Create residues that form a PPII helix
    residues = [
        create_mock_residue(phi=-75, psi=145),
        create_mock_residue(phi=-75, psi=145),
        create_mock_residue(phi=-75, psi=145),
    ]
    calculate_pp_helices(residues)

    # Check that the secondary structure and helix flags are set correctly
    assert residues[0].secondary_structure == StructureType.HELIX_PPII
    assert residues[0].helix_flags[HelixType.PP] == HelixPositionType.START
    assert residues[1].secondary_structure == StructureType.HELIX_PPII
    assert residues[1].helix_flags[HelixType.PP] == HelixPositionType.MIDDLE
    assert residues[2].secondary_structure == StructureType.HELIX_PPII
    assert residues[2].helix_flags[HelixType.PP] == HelixPositionType.END

def test_calculate_pp_helices_not_a_helix():
    """
    Tests that calculate_pp_helices does not assign a helix if the angles are wrong.
    """
    residues = [
        create_mock_residue(phi=0, psi=0),
        create_mock_residue(phi=0, psi=0),
        create_mock_residue(phi=0, psi=0),
    ]
    calculate_pp_helices(residues)
    assert residues[0].secondary_structure == StructureType.LOOP
    assert residues[1].secondary_structure == StructureType.LOOP
    assert residues[2].secondary_structure == StructureType.LOOP

def test_calculate_pp_helices_start_and_end():
    """
    Tests that calculate_pp_helices handles START_AND_END correctly.
    """
    residues = [
        create_mock_residue(phi=-75, psi=145),
        create_mock_residue(phi=-75, psi=145),
        create_mock_residue(phi=-75, psi=145),
    ]
    # Pre-set the end flag
    residues[0].helix_flags[HelixType.PP] = HelixPositionType.END
    calculate_pp_helices(residues)
    assert residues[0].helix_flags[HelixType.PP] == HelixPositionType.START_AND_END

def test_calculate_pp_helices_on_loop():
    """
    Tests that calculate_pp_helices correctly assigns a helix to residues
    that are initially in a LOOP state.
    """
    residues = [
        create_mock_residue(phi=-75, psi=145),
        create_mock_residue(phi=-75, psi=145),
        create_mock_residue(phi=-75, psi=145),
    ]
    residues[0].secondary_structure = StructureType.LOOP

    calculate_pp_helices(residues)

    assert residues[0].secondary_structure == StructureType.HELIX_PPII
