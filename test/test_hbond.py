import numpy as np
import pytest
import gzip
import re
from dsspy.io import read_cif
from dsspy.hbond import calculate_h_bonds


def parse_reference_dssp(filepath):
    """
    Parses a reference DSSP CIF file to extract H-bond information.
    Returns a dictionary where keys are residue numbers and values are
    dictionaries containing 'donor' and 'acceptor' H-bond lists.
    """
    hbonds = {}
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Find the _dssp_struct_bridge_pairs loop
    bridge_pairs_pattern = r'loop_\s+_dssp_struct_bridge_pairs\..*?\n(.*?)(?=\n#|\nloop_|\n_|\ndata_|\Z)'
    match = re.search(bridge_pairs_pattern, content, re.DOTALL)
    
    if not match:
        raise ValueError("Could not find _dssp_struct_bridge_pairs table in CIF file")
    
    # Extract the header to understand column positions
    header_pattern = r'loop_\s+((?:_dssp_struct_bridge_pairs\.\w+\s*\n)+)'
    header_match = re.search(header_pattern, content)
    
    if not header_match:
        raise ValueError("Could not parse _dssp_struct_bridge_pairs header")
    
    # Parse column names
    column_lines = header_match.group(1).strip().split('\n')
    columns = [line.strip().replace('_dssp_struct_bridge_pairs.', '') for line in column_lines]
    
    # Create column index mapping
    col_indices = {col: i for i, col in enumerate(columns)}
    
    # Required columns
    required_cols = ['label_seq_id', 'acceptor_1_label_seq_id', 'acceptor_1_energy', 
                    'acceptor_2_label_seq_id', 'acceptor_2_energy', 
                    'donor_1_label_seq_id', 'donor_1_energy', 
                    'donor_2_label_seq_id', 'donor_2_energy']
    
    for col in required_cols:
        if col not in col_indices:
            raise ValueError(f"Required column {col} not found in CIF file")
    
    # Parse the data lines
    data_section = match.group(1).strip()
    data_lines = [line.strip() for line in data_section.split('\n') if line.strip() and not line.strip().startswith('#')]
    
    for line in data_lines:
        if not line:
            continue
            
        # Split the line into fields, handling potential quoted values
        fields = []
        in_quotes = False
        current_field = ""
        
        i = 0
        while i < len(line):
            char = line[i]
            if char == '"' or char == "'":
                in_quotes = not in_quotes
                current_field += char
            elif char == ' ' and not in_quotes:
                if current_field:
                    fields.append(current_field)
                    current_field = ""
                # Skip multiple spaces
                while i + 1 < len(line) and line[i + 1] == ' ':
                    i += 1
            else:
                current_field += char
            i += 1
        
        if current_field:
            fields.append(current_field)
        
        if len(fields) < len(columns):
            continue
            
        # Extract residue number
        try:
            res_num = int(fields[col_indices['label_seq_id']])
        except (ValueError, IndexError):
            continue
            
        if res_num not in hbonds:
            hbonds[res_num] = {'donor': [], 'acceptor': []}
        
        # Process acceptor bonds (where this residue accepts H-bonds)
        for acceptor_col, energy_col in [('acceptor_1_label_seq_id', 'acceptor_1_energy'),
                                        ('acceptor_2_label_seq_id', 'acceptor_2_energy')]:
            try:
                partner_res_str = fields[col_indices[acceptor_col]]
                energy_str = fields[col_indices[energy_col]]
                
                if partner_res_str != '?' and energy_str != '?':
                    partner_res = int(partner_res_str)
                    energy = float(energy_str)
                    
                    if not np.isclose(energy, 0.0):
                        offset = partner_res - res_num
                        hbonds[res_num]['acceptor'].append({'offset': offset, 'energy': energy})
                        
            except (ValueError, IndexError):
                continue
        
        # Process donor bonds (where this residue donates H-bonds)
        for donor_col, energy_col in [('donor_1_label_seq_id', 'donor_1_energy'),
                                     ('donor_2_label_seq_id', 'donor_2_energy')]:
            try:
                partner_res_str = fields[col_indices[donor_col]]
                energy_str = fields[col_indices[energy_col]]
                
                if partner_res_str != '?' and energy_str != '?':
                    partner_res = int(partner_res_str)
                    energy = float(energy_str)
                    
                    if not np.isclose(energy, 0.0):
                        offset = partner_res - res_num
                        hbonds[res_num]['donor'].append({'offset': offset, 'energy': energy})
                        
            except (ValueError, IndexError):
                continue
    
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
