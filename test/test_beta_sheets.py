"""
Tests for the beta sheet calculation in the secondary_structure module.
"""
from unittest.mock import Mock
from dsspy.core import BridgeType, Residue, StructureType, BridgePartner
from dsspy.secondary_structure import calculate_beta_sheets

def create_mock_bio_residue(res_num, chain_id='A'):
    """Creates a mock biopython residue object."""
    mock_bio_res = Mock()
    mock_bio_res.get_id.return_value = (' ', res_num, ' ')
    mock_bio_res.get_full_id.return_value = ('', 0, chain_id, (' ', res_num, ' '))
    mock_bio_res.get_resname.return_value = 'ALA'
    # Make get_parent() return a mock object with an id attribute
    parent_mock = Mock()
    parent_mock.id = chain_id
    mock_bio_res.get_parent.return_value = parent_mock
    return mock_bio_res

def create_mock_residue(res_num, chain_id='A'):
    """Creates a mock residue object for testing."""
    bio_res = create_mock_bio_residue(res_num, chain_id)
    res = Residue(bio_res, res_num)
    res.hbond_acceptor = []
    res.hbond_donor = []
    res.beta_partner = [BridgePartner(None, 0, False), BridgePartner(None, 0, False)]
    res.sheet = 0
    res.secondary_structure = StructureType.LOOP
    return res

def test_calculate_beta_sheets_empty_list():
    """
    Tests that calculate_beta_sheets handles an empty list of residues.
    """
    calculate_beta_sheets([])
    # No assertion needed, just testing that it doesn't crash.

def test_antiparallel_ladder_merging():
    """
    Tests that two antiparallel ladders that are close together are merged.
    """
    res = {i: create_mock_residue(i) for i in range(1, 12)}
    for i in range(1, 11):
        res[i].next_residue = res[i+1]
        res[i+1].prev_residue = res[i]
    
    residues = list(res.values())

    # Bridge 1: 2 <> 9 (antiparallel). Pattern: (c,d) and (f,a)
    # _test_bridge(res2, res9): a=1, b=2, c=3, d=8, e=9, f=10
    # H-bonds: res3 accepts from res8, res10 accepts from res1
    res[3].hbond_acceptor.append(Mock(residue=res[8], energy=-2.0))
    res[10].hbond_acceptor.append(Mock(residue=res[1], energy=-2.0))

    # Bridge 2: 4 <> 7 (antiparallel). Pattern: (e,b) and (b,e)
    # _test_bridge(res4, res7): a=3, b=4, c=5, d=6, e=7, f=8
    # H-bonds: res7 accepts from res4, res4 accepts from res7
    res[7].hbond_acceptor.append(Mock(residue=res[4], energy=-2.0))
    res[4].hbond_acceptor.append(Mock(residue=res[7], energy=-2.0))

    calculate_beta_sheets(residues)

    # Check that the bridges were merged into a single sheet
    assert res[2].sheet == 1
    assert res[4].sheet == 1
    assert res[7].sheet == 1
    assert res[9].sheet == 1
    assert res[2].beta_partner[0].parallel is False
    assert res[4].beta_partner[0].parallel is False

def test_residue_in_multiple_sheets():
    """
    Tests that a residue can be part of multiple beta sheets.
    """
    res = {i: create_mock_residue(i) for i in range(1, 12)}
    for i in range(1, 11):
        res[i].next_residue = res[i+1]
        res[i+1].prev_residue = res[i]
    
    residues = list(res.values())

    # Bridge 1: 2 <> 5 (antiparallel)
    res[5].hbond_acceptor.append(Mock(residue=res[2], energy=-2.0))
    res[2].hbond_acceptor.append(Mock(residue=res[5], energy=-2.0))
    
    # Bridge 2: 8 <> 10 (antiparallel) - separate from bridge 1
    res[10].hbond_acceptor.append(Mock(residue=res[8], energy=-2.0))
    res[8].hbond_acceptor.append(Mock(residue=res[10], energy=-2.0))

    calculate_beta_sheets(residues)

    # Now, add a third bridge to link the two sheets
    # Bridge 3: 5 <> 8
    res[8].hbond_acceptor.append(Mock(residue=res[5], energy=-2.0))
    res[5].hbond_acceptor.append(Mock(residue=res[8], energy=-2.0))

    # Recalculate with the new bridge, which should merge the sheets
    calculate_beta_sheets(residues)

    # Now res5 should have two partners
    assert res[5].beta_partner[0].residue is not None
    assert res[5].beta_partner[1].residue is not None
    partners = {p.residue.number for p in res[5].beta_partner if p.residue}
    assert partners == {2, 8}

def test_residue_already_in_strand():
    """
    Tests that if a residue is already part of a STRAND, it is not demoted
    to a BETA_BRIDGE.
    """
    res = {i: create_mock_residue(i) for i in range(1, 11)}
    for i in range(1, 10):
        res[i].next_residue = res[i+1]
        res[i+1].prev_residue = res[i]
    
    residues = list(res.values())

    # Create a strand of length 2 (res2, res3)
    # Bridge 1: 2 <> 7 (antiparallel)
    res[2].hbond_acceptor.append(Mock(residue=res[7], energy=-2.0))
    res[7].hbond_acceptor.append(Mock(residue=res[2], energy=-2.0))
    # Bridge 2: 3 <> 6 (antiparallel)
    res[3].hbond_acceptor.append(Mock(residue=res[6], energy=-2.0))
    res[6].hbond_acceptor.append(Mock(residue=res[3], energy=-2.0))
    
    calculate_beta_sheets(residues)

    # res[2] and res[3] should now be in a STRAND
    assert res[2].secondary_structure == StructureType.STRAND
    assert res[3].secondary_structure == StructureType.STRAND

    # Now, add a new, single-residue bridge that involves residue 3
    # Bridge 3: 3 <> 9 (antiparallel)
    res[3].hbond_acceptor.append(Mock(residue=res[9], energy=-2.0))
    res[9].hbond_acceptor.append(Mock(residue=res[3], energy=-2.0))

    calculate_beta_sheets(residues)

    # Check that residue 3 is still a STRAND
    assert res[3].secondary_structure == StructureType.STRAND
