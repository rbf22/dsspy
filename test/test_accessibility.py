"""
Tests for accessibility calculation.
"""

import re
import pytest
from dsspy.io import read_cif
from dsspy.accessibility import calculate_accessibility


def parse_reference_accessibility(filepath):
    """
    Parses a reference DSSP file in mmCIF format to extract accessibility values.
    Returns a dictionary where keys are residue numbers and values are the accessibility.
    """
    accessibilities = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    in_summary_loop = False
    seq_id_index = -1
    accessibility_index = -1
    header_count = 0

    for line in lines:
        line = line.strip()
        if line.startswith('loop_'):
            in_summary_loop = False # Reset on new loop
            header_count = 0
        elif line.startswith('_dssp_struct_summary.'):
            in_summary_loop = True
            if 'label_seq_id' in line:
                seq_id_index = header_count
            elif 'accessibility' in line:
                accessibility_index = header_count
            header_count += 1
        elif in_summary_loop and not line.startswith('#'):
            if line.startswith('loop_'): # Should not happen with the new file
                in_summary_loop = False
                continue

            parts = re.split(r'\s+', line)
            if len(parts) > max(seq_id_index, accessibility_index):
                res_num_str = parts[seq_id_index]
                acc_str = parts[accessibility_index]
                if acc_str not in {".", "?"}:
                    accessibilities[int(res_num_str)] = float(acc_str)
        elif not line.startswith("_"):
            header_count = 0
    return accessibilities

def test_calculate_accessibility_comparative():
    """
    Tests the calculate_accessibility function by comparing its output to a
    reference DSSP file.
    """
    # 1. Run dsspy's accessibility calculation
    with open('test/reference_data/1cbs.cif', 'rt', encoding='utf-8') as f:
        residues, _ = read_cif(f)
    calculate_accessibility(residues)

    # 2. Parse the reference DSSP file
    reference_accessibilities = parse_reference_accessibility(
        'test/reference_data/1cbs-dssp.cif'
    )

    # 3. Compare the results
    for res in residues:
        if res.number in reference_accessibilities:
            ref_acc = reference_accessibilities[res.number]
            # It seems the accessibility is an integer in the reference file,
            # but a float in our calculation.
            # Let's round our result for comparison.
            calculated_acc = round(res.accessibility)
            assert calculated_acc == pytest.approx(ref_acc, abs=1)
