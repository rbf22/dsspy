import numpy as np
import pytest
import gzip
import re
from dsspy.io import read_cif
from dsspy.hbond import calculate_h_bonds


def parse_reference_dssp(filepath):
    """
    Parses a reference DSSP file to extract H-bond information.
    Returns a dictionary where keys are residue numbers and values are
    dictionaries containing 'donor' and 'acceptor' H-bond lists.
    """
    hbonds = {}
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Find the start of the residue data
    data_start_index = -1
    for i, line in enumerate(lines):
        if line.strip().startswith('#  RESIDUE AA STRUCTURE'):
            data_start_index = i + 1
            break

    if data_start_index == -1:
        raise ValueError("Could not find data section in reference DSSP file")

    for line in lines[data_start_index:]:
        if len(line) < 100: # Skip malformed or summary lines
            continue

        res_num = int(line[0:5].strip())
        hbonds[res_num] = {'donor': [], 'acceptor': []}

        # Extract H-bonds
        # N-H-->O (donor) columns are at indices 38 and 50
        # O-->H-N (acceptor) columns are at indices 62 and 74
        donor1_str = line[38:49].strip()
        donor2_str = line[50:61].strip()
        acceptor1_str = line[62:73].strip()
        acceptor2_str = line[74:85].strip()

        for d_str in [donor1_str, donor2_str]:
            if d_str:
                offset, energy_str = d_str.split(',')
                parts = re.findall(r"[-+]?\d*\.\d+", energy_str)
                if parts:
                    energy = float(parts[0])
                else:
                    continue
                if not np.isclose(energy, 0.0):
                    hbonds[res_num]['donor'].append({'offset': int(offset), 'energy': energy})

        for a_str in [acceptor1_str, acceptor2_str]:
            if a_str:
                offset, energy_str = a_str.split(',')
                parts = re.findall(r"[-+]?\d*\.\d+", energy_str)
                if parts:
                    energy = float(parts[0])
                else:
                    continue
                if not np.isclose(energy, 0.0):
                    hbonds[res_num]['acceptor'].append({'offset': int(offset), 'energy': energy})

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
    reference_hbonds = parse_reference_dssp('test/reference_data/1cbs.dssp')

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
