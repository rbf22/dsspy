"""
Utility functions for the dsspy package.
"""

def get_sheet_label(sheet_number):
    """
    Converts a sheet number to a sheet label (A-Z, a-z).
    """
    if sheet_number == 0:
        return ' '
    if 1 <= sheet_number <= 26:
        return chr(ord('A') + sheet_number - 1)
    if 27 <= sheet_number <= 52:
        return chr(ord('a') + sheet_number - 27)
    return '?'
