"""
Tests for the core data structures.
"""

from unittest.mock import Mock
from dsspy.core import HBond, Residue, BridgePartner, Bridge, BridgeType

def test_hbond_repr():
    """Tests the __repr__ of the HBond class."""
    mock_residue = Mock()
    mock_residue.id = "A/1"
    hbond_with_residue = HBond(mock_residue, -0.5)
    assert repr(hbond_with_residue) == "<HBond to A/1 with energy -0.5>"

    hbond_without_residue = HBond(None, 0)
    assert repr(hbond_without_residue) == "<HBond to None with energy 0>"

def test_residue_repr():
    """Tests the __repr__ of the Residue class."""
    mock_bio_res = Mock()
    mock_bio_res.get_resname.return_value = 'ALA'
    mock_bio_res.get_id.return_value = (' ', 1, ' ')
    res = Residue(mock_bio_res, 1)
    assert repr(res) == "<Residue ALA id=(' ', 1, ' ')>"

def test_bridge_partner_init():
    """Tests the __init__ of the BridgePartner class."""
    mock_residue = Mock()
    bp = BridgePartner(mock_residue, ladder=1, parallel=True)
    assert bp.residue == mock_residue
    assert bp.ladder == 1
    assert bp.parallel is True

def test_bridge_partner_repr():
    """Tests the __repr__ of the BridgePartner class."""
    mock_residue = Mock()
    mock_residue.id = "A/1"
    bp = BridgePartner(mock_residue)
    assert repr(bp) == "<BridgePartner to A/1>"

def test_bridge_repr():
    """Tests the __repr__ of the Bridge class."""
    mock_res1 = Mock()
    mock_res1.id = "A/1"
    mock_res2 = Mock()
    mock_res2.id = "A/2"
    bridge = Bridge(mock_res1, mock_res2, BridgeType.PARALLEL)
    assert repr(bridge) == "<Bridge A/1-A/2 type=BridgeType.PARALLEL>"
