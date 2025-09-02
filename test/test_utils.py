"""
Tests for the utils module.
"""

from dsspy.utils import get_sheet_label

def test_get_sheet_label():
    """
    Tests the get_sheet_label function.
    """
    assert get_sheet_label(0) == ' '
    assert get_sheet_label(1) == 'A'
    assert get_sheet_label(26) == 'Z'
    assert get_sheet_label(27) == 'a'
    assert get_sheet_label(52) == 'z'
    assert get_sheet_label(53) == '?'
    assert get_sheet_label(-1) == '?'
