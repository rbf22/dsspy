import numpy as np
import pytest
import gzip
import re
from dsspy.io import read_cif
from dsspy.hbond import calculate_h_bonds


def parse_reference_dssp(filepath):
    hbonds = {}
    in_loop = False
    columns = []
    col_indices = {}

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Detect start of the bridge_pairs loop
            if line.startswith("loop_"):
                in_loop = False  # reset
                columns = []
                continue

            if line.startswith("_dssp_struct_bridge_pairs."):
                in_loop = True
                columns.append(line.replace("_dssp_struct_bridge_pairs.", ""))
                continue

            # If we are inside the loop and see a new header or section, stop
            if in_loop and (line.startswith("_") or line.startswith("loop_") or line.startswith("data_") or line.startswith("#")):
                in_loop = False
                continue

            # Parse data rows
            if in_loop:
                if not col_indices:
                    col_indices = {col: i for i, col in enumerate(columns)}

                fields = line.split()
                if len(fields) < len(columns):
                    continue

                try:
                    res_num = int(fields[col_indices['label_seq_id']])
                except Exception:
                    continue

                if res_num not in hbonds:
                    hbonds[res_num] = {'donor': [], 'acceptor': []}

                # acceptors
                for acc, ene in [
                    ('acceptor_1_label_seq_id', 'acceptor_1_energy'),
                    ('acceptor_2_label_seq_id', 'acceptor_2_energy'),
                ]:
                    try:
                        partner = fields[col_indices[acc]]
                        energy = fields[col_indices[ene]]
                        if partner != '?' and energy != '?':
                            offset = int(partner) - res_num
                            hbonds[res_num]['acceptor'].append({
                                'offset': offset,
                                'energy': float(energy)
                            })
                    except Exception:
                        pass

                # donors
                for don, ene in [
                    ('donor_1_label_seq_id', 'donor_1_energy'),
                    ('donor_2_label_seq_id', 'donor_2_energy'),
                ]:
                    try:
                        partner = fields[col_indices[don]]
                        energy = fields[col_indices[ene]]
                        if partner != '?' and energy != '?':
                            offset = int(partner) - res_num
                            hbonds[res_num]['donor'].append({
                                'offset': offset,
                                'energy': float(energy)
                            })
                    except Exception:
                        pass

    return hbonds

def test_calculate_h_bonds_comparative():
    """
    Tests the calculate_h_bonds function by comparing its output to a
    reference DSSP file.
    """
    # 1. Run dsspy's H-bond calculation
    with gzip.open('test/reference_data/1cbs.cif.gz', 'rt') as f:
        residues, _ = read_cif(f)
    calculate_h_bonds(residues)

    # Create a map of residue number to residue object for easy lookup
    res_map = {res.id: res for res in residues}

    # 2. Parse the reference DSSP file
    reference_hbonds = parse_reference_dssp('test/reference_data/1cbs-dssp.cif')

    # 3. Compare the results
    assert len(residues) == len(reference_hbonds)

    for i, res in enumerate(residues, 1):
        # Make sure our enumeration index matches the residue's internal number
        assert i == res.number

        ref_bonds = reference_hbonds.get(i)
        if not ref_bonds:
            continue

        # Compare donor bonds (N-H --> O)
        # NOTE: C++ code uses counter-intuitive list names.
        # hbond_acceptor is for bonds where `res` is the DONOR.
        dsspy_donors = sorted([
            {'offset': h.residue.number - res.number, 'energy': h.energy}
            for h in res.hbond_acceptor if h.residue is not None
        ], key=lambda x: x['offset'])

        ref_donors = sorted(ref_bonds['donor'], key=lambda x: x['offset'])

        if len(dsspy_donors) != len(ref_donors):
            print(f"Residue {i}: Mismatch in number of donor bonds")
            print(f"  dsspy: {dsspy_donors}")
            print(f"  ref:   {ref_donors}")
        assert len(dsspy_donors) == len(ref_donors)

        for dsspy_d, ref_d in zip(dsspy_donors, ref_donors):
            if dsspy_d['offset'] != ref_d['offset']:
                print(f"Residue {i}: Mismatch in donor bond offset")
                print(f"  dsspy: {dsspy_d}")
                print(f"  ref:   {ref_d}")
            assert dsspy_d['offset'] == ref_d['offset']
            assert dsspy_d['energy'] == pytest.approx(ref_d['energy'], abs=1e-3)

        # Compare acceptor bonds (O --> H-N)
        # NOTE: C++ code uses counter-intuitive list names.
        # hbond_donor is for bonds where `res` is the ACCEPTOR.
        dsspy_acceptors = sorted([
            {'offset': h.residue.number - res.number, 'energy': h.energy}
            for h in res.hbond_donor if h.residue is not None
        ], key=lambda x: x['offset'])

        ref_acceptors = sorted(ref_bonds['acceptor'], key=lambda x: x['offset'])

        if len(dsspy_acceptors) != len(ref_acceptors):
            print(f"Residue {i}: Mismatch in number of acceptor bonds")
            print(f"  dsspy: {dsspy_acceptors}")
            print(f"  ref:   {ref_acceptors}")
        assert len(dsspy_acceptors) == len(ref_acceptors)

        for dsspy_a, ref_a in zip(dsspy_acceptors, ref_acceptors):
            if dsspy_a['offset'] != ref_a['offset']:
                print(f"Residue {i}: Mismatch in acceptor bond offset")
                print(f"  dsspy: {dsspy_a}")
                print(f"  ref:   {ref_a}")
            assert dsspy_a['offset'] == ref_a['offset']
            assert dsspy_a['energy'] == pytest.approx(ref_a['energy'], abs=1e-3)
