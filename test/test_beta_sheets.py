"""
Tests for the beta sheet calculation in the secondary_structure module.
"""
from unittest.mock import Mock
from dsspy.core import BridgeType, Residue, StructureType, BridgePartner
from dsspy.secondary_structure import calculate_beta_sheets

def create_mock_bio_residue(res_num, chain_id='A', parent_mock=None):
    """Creates a mock biopython residue object."""
    mock_bio_res = Mock()
    mock_bio_res.get_id.return_value = (' ', res_num, ' ')
    mock_bio_res.get_full_id.return_value = ('', 0, chain_id, (' ', res_num, ' '))
    mock_bio_res.get_resname.return_value = 'ALA'
    # Make get_parent() return a mock object with an id attribute
    if parent_mock is None:
        parent_mock = Mock()
        parent_mock.id = chain_id
    mock_bio_res.get_parent.return_value = parent_mock
    return mock_bio_res

def create_mock_residue(res_num, chain_id='A', parent_mock=None):
    """Creates a mock residue object for testing."""
    bio_res = create_mock_bio_residue(res_num, chain_id, parent_mock)
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
    parent_mock = Mock()
    parent_mock.id = 'A'
    res = {i: create_mock_residue(i, parent_mock=parent_mock) for i in range(1, 12)}
    for i in range(1, 11):
        res[i].next_residue = res[i+1]
        res[i+1].prev_residue = res[i]
    
    residues = list(res.values())

    # Bridge 1: 2 <> 9 (antiparallel). Pattern: (c,d) and (f,a)
    # _test_bridge(res2, res9): a=1, b=2, c=3, d=8, e=9, f=10
    # H-bonds: res3 accepts from res8, res10 accepts from res1
    # This means 8 -> 3 and 1 -> 10.
    # The donor's hbond_acceptor list contains the acceptor.
    res[8].hbond_acceptor.append(Mock(residue=res[3], energy=-2.0))
    res[1].hbond_acceptor.append(Mock(residue=res[10], energy=-2.0))

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
    parent_mock = Mock()
    parent_mock.id = 'A'
    res = {i: create_mock_residue(i, parent_mock=parent_mock) for i in range(1, 12)}
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
    parent_mock = Mock()
    parent_mock.id = 'A'
    res = {i: create_mock_residue(i, parent_mock=parent_mock) for i in range(1, 11)}
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

def test_parallel_bridge():
    """
    Tests that parallel bridges are correctly identified and merged.
    """
    parent_mock = Mock()
    parent_mock.id = 'A'
    res = {i: create_mock_residue(i, parent_mock=parent_mock) for i in range(1, 10)}
    for i in range(1, 9):
        res[i].next_residue = res[i+1]
        res[i+1].prev_residue = res[i]

    residues = list(res.values())

    # Bridge 1: 2 <> 5
    # Pattern: (i-1, j) and (j, i+1) => (1,5) and (5,3)
    res[1].hbond_acceptor.append(Mock(residue=res[5], energy=-2.0))
    res[5].hbond_acceptor.append(Mock(residue=res[3], energy=-2.0))

    # Bridge 2: 3 <> 6
    # Pattern: (i-1, j) and (j, i+1) => (2,6) and (6,4)
    res[2].hbond_acceptor.append(Mock(residue=res[6], energy=-2.0))
    res[6].hbond_acceptor.append(Mock(residue=res[4], energy=-2.0))

    calculate_beta_sheets(residues)

    # Check that sheets are formed and merged
    assert res[2].sheet == 1
    assert res[3].sheet == 1
    assert res[5].sheet == 1
    assert res[6].sheet == 1

    # Check partners
    assert res[2].beta_partner[0].residue == res[5]
    assert res[3].beta_partner[0].residue == res[6]
    assert res[5].beta_partner[0].residue == res[2]
    assert res[6].beta_partner[0].residue == res[3]
    assert res[2].beta_partner[0].parallel is True

def test_bridge_across_chains():
    """
    Tests that a bridge is not formed between residues in different chains.
    """
    # Create residues in two different chains
    parent_A = Mock()
    parent_A.id = 'A'
    res_A = {i: create_mock_residue(i, chain_id='A', parent_mock=parent_A) for i in range(1, 4)}

    parent_B = Mock()
    parent_B.id = 'B'
    res_B = {i: create_mock_residue(i, chain_id='B', parent_mock=parent_B) for i in range(4, 7)}

    # Link residues within each chain
    for i in range(1, 3):
        res_A[i].next_residue = res_A[i+1]
        res_A[i+1].prev_residue = res_A[i]
    for i in range(4, 6):
        res_B[i].next_residue = res_B[i+1]
        res_B[i+1].prev_residue = res_B[i]

    residues = list(res_A.values()) + list(res_B.values())

    # Try to create a bridge between res 2 (chain A) and 5 (chain B)
    res_A[1].hbond_acceptor.append(Mock(residue=res_B[5], energy=-2.0))
    res_B[5].hbond_acceptor.append(Mock(residue=res_A[3], energy=-2.0))

    calculate_beta_sheets(residues)

    # Assert that no sheet is formed
    assert res_A[2].sheet == 0
    assert res_B[5].sheet == 0

def test_bridge_with_chain_gap():
    """
    Tests that a bridge is not formed if there is a gap in the chain.
    """
    parent_mock = Mock()
    parent_mock.id = 'A'
    res = {i: create_mock_residue(i, parent_mock=parent_mock) for i in range(1, 7)}

    # Create a gap between residue 2 and 3
    res[1].next_residue = res[2]
    res[2].prev_residue = res[1]
    # No link from 2 to 3
    res[3].next_residue = res[4]
    res[4].prev_residue = res[3]
    res[4].next_residue = res[5]
    res[5].prev_residue = res[4]
    res[5].next_residue = res[6]
    res[6].prev_residue = res[5]

    residues = list(res.values())

    # Try to create a bridge between res 2 and 5, which would require
    # a continuous chain from 1-3 and 4-6. The 1-3 chain is broken.
    res[1].hbond_acceptor.append(Mock(residue=res[5], energy=-2.0))
    res[5].hbond_acceptor.append(Mock(residue=res[3], energy=-2.0))

    calculate_beta_sheets(residues)

    # Assert that no sheet is formed because of the break
    assert res[2].sheet == 0
    assert res[5].sheet == 0
