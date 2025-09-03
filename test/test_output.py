"""
Tests for the output module.
"""

from unittest.mock import Mock
from dsspy.output import _format_header, format_dssp_line
from dsspy.io import read_cif
from dsspy.hbond import calculate_h_bonds
from dsspy.secondary_structure import calculate_beta_sheets
from dsspy.utils import get_sheet_label
from dsspy.core import HelixType, HelixPositionType, StructureType, Residue

def test_format_header():
    """
    Tests the _format_header function with a sample header dictionary.
    """
    header = {
        'head': 'test protein',
        'deposition_date': '1999-10-26',
        'idcode': 'XXXX',
        'compound': [{'mol_id': '1', 'molecule': 'test protein', 'chain': 'A'}],
        'source': [{'mol_id': '1', 'organism_scientific': 'homo sapiens'}],
        'author': 'J. Doe'
    }

    header_str = _format_header(header)

    assert "HEADER    TEST PROTEIN" in header_str
    assert "1999-10-26" in header_str
    assert "XXXX" in header_str
    assert "COMPND    MOL_ID: 1; MOLECULE: test protein; CHAIN: A" in header_str
    assert "SOURCE    MOL_ID: 1; ORGANISM_SCIENTIFIC: homo sapiens" in header_str
    assert "AUTHOR    J. Doe" in header_str

def test_format_header_empty():
    """
    Tests the _format_header function with an empty header dictionary.
    """
    header_str = _format_header({})
    assert "HEADER" in header_str
    assert "COMPND    " in header_str
    assert "SOURCE    " in header_str
    assert "AUTHOR" in header_str

    # Test with empty compound and source lists
    header_str = _format_header({'compound': [], 'source': []})
    assert "COMPND    " in header_str
    assert "SOURCE    " in header_str

    # Test with non-empty dictionaries in the lists
    header_str = _format_header({'compound': [{'mol_id': '1'}], 'source': [{'mol_id': '1'}]})
    assert "COMPND    MOL_ID: 1" in header_str
    assert "SOURCE    MOL_ID: 1" in header_str


def test_format_dssp_line_with_sheets():
    """
    Tests that the sheet and bridge labels are correctly formatted in the output line.
    This is an integration test for calculate_beta_sheets and format_dssp_line.
    """
    with open('test/reference_data/1cbs.cif', 'rt', encoding='utf-8') as f:
        residues, _ = read_cif(f)

    # Run the algorithms that calculate sheet information
    calculate_h_bonds(residues)
    calculate_beta_sheets(residues)

    # Find a residue that is part of a bridge
    bridge_res = None
    for res in residues:
        if res.beta_partner[0].residue is not None:
            bridge_res = res
            break

    assert bridge_res is not None, "No residue with a bridge partner found"
    assert bridge_res.sheet > 0

    line = format_dssp_line(bridge_res)

    # Check that the sheet and bridge labels are formatted correctly
    bp1 = bridge_res.beta_partner[0].residue.number
    bp2 = bridge_res.beta_partner[1].residue.number if bridge_res.beta_partner[1].residue else 0

    sheet_label = get_sheet_label(bridge_res.sheet)
    ladder_label = chr(ord('a') + bridge_res.beta_partner[0].ladder % 26)

    expected_str = f"{bp1:>4d}{bp2:>4d}{ladder_label}{sheet_label}"

    assert expected_str in line

def create_mock_residue():
    """Creates a mock residue object for testing."""
    mock_bio_res = Mock()
    mock_bio_res.get_id.return_value = (' ', 1, ' ')
    mock_bio_res.get_full_id.return_value = ('', 0, 'A', (' ', 1, ' '))
    mock_bio_res.get_resname.return_value = 'ALA'

    # Mock the __getitem__ method to allow dictionary-style access for atoms
    mock_atom = Mock()
    mock_atom.get_coord.return_value = (0.0, 0.0, 0.0)
    mock_bio_res.__getitem__ = Mock(return_value=mock_atom)

    # Create a mock Residue object
    # We can't just mock the whole thing because format_dssp_line expects a real Residue object
    # so it can access attributes directly.
    res = Residue(mock_bio_res, 1)
    res.secondary_structure = StructureType.LOOP
    res.hbond_acceptor = [Mock(), Mock()]
    res.hbond_acceptor[0].residue = None
    res.hbond_acceptor[1].residue = None
    res.hbond_donor = [Mock(), Mock()]
    res.hbond_donor[0].residue = None
    res.hbond_donor[1].residue = None
    res.beta_partner = [Mock(), Mock()]
    res.beta_partner[0].residue = None
    res.beta_partner[1].residue = None
    res.sheet = 0
    res.bend = False
    res.alpha = 0.0
    res.accessibility = 0.0
    res.tco = 0.0
    res.kappa = 0.0
    res.phi = 0.0
    res.psi = 0.0

    return res

def test_format_dssp_line_helix_flags():
    """
    Tests that the helix flags are correctly formatted in the output line.
    """
    res = create_mock_residue()

    # Test HelixPositionType.START
    res.helix_flags[HelixType.THREE_TEN] = HelixPositionType.START
    line = format_dssp_line(res)
    assert " >   " in line
    res.helix_flags[HelixType.THREE_TEN] = HelixPositionType.NONE

    # Test HelixPositionType.END
    res.helix_flags[HelixType.ALPHA] = HelixPositionType.END
    line = format_dssp_line(res)
    assert "   < " in line
    res.helix_flags[HelixType.ALPHA] = HelixPositionType.NONE

    # Test HelixPositionType.START_AND_END
    res.helix_flags[HelixType.ALPHA] = HelixPositionType.START_AND_END
    line = format_dssp_line(res)
    assert "   X " in line
    res.helix_flags[HelixType.ALPHA] = HelixPositionType.NONE

    # Test HelixPositionType.MIDDLE for PI helix
    res.helix_flags[HelixType.PI] = HelixPositionType.MIDDLE
    line = format_dssp_line(res)
    assert "  5  " in line
    res.helix_flags[HelixType.PI] = HelixPositionType.NONE

    # Test HelixPositionType.MIDDLE for PP helix
    res.helix_flags[HelixType.PP] = HelixPositionType.MIDDLE
    line = format_dssp_line(res)
    assert "   P " in line
    res.helix_flags[HelixType.PP] = HelixPositionType.NONE

def test_format_header_with_compound_and_source():
    """
    Tests the _format_header function with compound and source information.
    """
    header = {
        'compound': [{'mol_id': '1', 'molecule': 'test protein'}],
        'source': [{'mol_id': '1', 'organism_scientific': 'homo sapiens'}]
    }
    header_str = _format_header(header)
    assert "COMPND    MOL_ID: 1; MOLECULE: test protein" in header_str
    assert "SOURCE    MOL_ID: 1; ORGANISM_SCIENTIFIC: homo sapiens" in header_str
