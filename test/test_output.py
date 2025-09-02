from dsspy.output import _format_header, format_dssp_line
from dsspy.io import read_cif
from dsspy.hbond import calculate_h_bonds
from dsspy.secondary_structure import calculate_beta_sheets

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

def test_format_dssp_line_with_sheets():
    """
    Tests that the sheet and bridge labels are correctly formatted in the output line.
    This is an integration test for calculate_beta_sheets and format_dssp_line.
    """
    with open('test/reference_data/1cbs.cif', 'rt') as f:
        residues, structure = read_cif(f)

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

    # Helper to get sheet label (A-Z, a-z)
    def _get_sheet_label(sheet_number):
        if sheet_number == 0:
            return ' '
        if 1 <= sheet_number <= 26:
            return chr(ord('A') + sheet_number - 1)
        if 27 <= sheet_number <= 52:
            return chr(ord('a') + sheet_number - 27)
        return '?'

    sheet_label = _get_sheet_label(bridge_res.sheet)
    ladder_label = chr(ord('a') + bridge_res.beta_partner[0].ladder % 26)

    expected_str = f"{bp1:>4d}{bp2:>4d}{ladder_label}{sheet_label}"

    assert expected_str in line
