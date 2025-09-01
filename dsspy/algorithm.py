import numpy as np
from dsspy.core import Residue

K_MIN_HBOND_ENERGY = -9.9
K_MAX_HBOND_ENERGY = -0.5
K_COUPLING_CONSTANT = -27.888

def distance(p1, p2):
    return np.linalg.norm(p1 - p2)

def calculate_hbond_energy(donor, acceptor):
    """
    Calculates the H-bond energy between a donor and an acceptor residue.
    """
    if donor.resname == 'PRO':
        return 0.0

    dist_ho = distance(donor.h_coord, acceptor.biopython_residue['O'].get_coord())
    dist_hc = distance(donor.h_coord, acceptor.biopython_residue['C'].get_coord())
    dist_nc = distance(donor.biopython_residue['N'].get_coord(), acceptor.biopython_residue['C'].get_coord())
    dist_no = distance(donor.biopython_residue['N'].get_coord(), acceptor.biopython_residue['O'].get_coord())

    energy = K_COUPLING_CONSTANT / dist_ho - K_COUPLING_CONSTANT / dist_hc + \
             K_COUPLING_CONSTANT / dist_nc - K_COUPLING_CONSTANT / dist_no

    energy = round(energy * 1000) / 1000

    if energy < K_MIN_HBOND_ENERGY:
        energy = K_MIN_HBOND_ENERGY

    return energy

def calculate_hbonds(residues):
    """
    Calculates all H-bonds for a list of residues.
    """
    # First, assign hydrogen positions
    for res in residues:
        res.assign_hydrogen()

    # This is a naive N^2 loop. The C++ code has an optimization to only
    # check pairs of residues that are close to each other.
    for i, r1 in enumerate(residues):
        for j, r2 in enumerate(residues):
            if i == j:
                continue

            # H-bond from r1 to r2
            energy1 = calculate_hbond_energy(r1, r2)
            if energy1 < K_MAX_HBOND_ENERGY:
                r1.hbond_acceptor.append({'residue': r2, 'energy': energy1})
                r2.hbond_donor.append({'residue': r1, 'energy': energy1})

            # H-bond from r2 to r1
            energy2 = calculate_hbond_energy(r2, r1)
            if energy2 < K_MAX_HBOND_ENERGY:
                r2.hbond_acceptor.append({'residue': r1, 'energy': energy2})
                r1.hbond_donor.append({'residue': r2, 'energy': energy2})

    # Sort the h-bonds by energy
    for res in residues:
        res.hbond_acceptor.sort(key=lambda x: x['energy'])
        res.hbond_donor.sort(key=lambda x: x['energy'])


def test_bond(r1, r2):
    """
    Checks if there is an H-bond from r1 to r2.
    """
    for acceptor in r1.hbond_acceptor:
        if acceptor['residue'] == r2 and acceptor['energy'] < K_MAX_HBOND_ENERGY:
            return True
    return False

def test_bridge(r1, r2):
    """
    Tests if two residues form a parallel or anti-parallel bridge.
    """
    # Based on the C++ code's TestBridge function
    a = r1.prev_residue
    b = r1
    c = r1.next_residue
    d = r2.prev_residue
    e = r2
    f = r2.next_residue

    if not (a and c and d and f):
        return None

    # Parallel bridge
    if (test_bond(c, e) and test_bond(e, a)) or \
       (test_bond(f, b) and test_bond(b, d)):
        return 'parallel'

    # Anti-parallel bridge
    if (test_bond(c, d) and test_bond(f, a)) or \
       (test_bond(e, b) and test_bond(b, e)):
        return 'antiparallel'

    return None

from dsspy.core import Bridge, StructureType

def calculate_beta_sheets(residues):
    """
    Calculates beta sheets from the H-bond information.
    """
    # 1. Find all bridges
    bridges = []
    for i, r1 in enumerate(residues):
        for j, r2 in enumerate(residues):
            # Residues in a bridge must be separated by at least 2 residues
            if abs(i - j) < 2:
                continue

            bridge_type = test_bridge(r1, r2)
            if bridge_type:
                bridges.append(Bridge(r1, r2, bridge_type))

    # 2. Group bridges into ladders
    ladders = []
    bridges.sort(key=lambda b: (b.r1.number, b.r2.number))

    for bridge in bridges:
        found_ladder = False
        for ladder in ladders:
            last_bridge = ladder[-1]
            if bridge.type == last_bridge.type and \
               bridge.r1.number == last_bridge.r1.number + 1:
                if bridge.type == 'parallel' and bridge.r2.number == last_bridge.r2.number + 1:
                    ladder.append(bridge)
                    found_ladder = True
                    break
                elif bridge.type == 'antiparallel' and bridge.r2.number == last_bridge.r2.number - 1:
                    ladder.append(bridge)
                    found_ladder = True
                    break
        if not found_ladder:
            ladders.append([bridge])

    # 3. Group ladders into sheets
    sheets = []
    ladders_in_sheets = set()
    for i, ladder1 in enumerate(ladders):
        if i in ladders_in_sheets:
            continue

        new_sheet = {i}
        q = [i]
        head = 0
        while head < len(q):
            current_ladder_idx = q[head]
            head += 1

            for j, ladder2 in enumerate(ladders):
                if j in new_sheet:
                    continue

                # Check if ladders are linked
                linked = False
                for b1 in ladders[current_ladder_idx]:
                    for b2 in ladder2:
                        if b1.r1 == b2.r1 or b1.r1 == b2.r2 or \
                           b1.r2 == b2.r1 or b1.r2 == b2.r2:
                            linked = True
                            break
                    if linked:
                        break

                if linked:
                    new_sheet.add(j)
                    q.append(j)

        sheets.append([ladders[k] for k in new_sheet])
        ladders_in_sheets.update(new_sheet)

    # 4. Assign sheet and structure information to residues
    sheet_nr = 1
    for sheet in sheets:
        for ladder in sheet:
            for bridge in ladder:
                for res in [bridge.r1, bridge.r2]:
                    if res.secondary_structure != StructureType.STRAND:
                        res.secondary_structure = StructureType.BETA_BRIDGE
                    if len(ladder) > 1:
                        res.secondary_structure = StructureType.STRAND
                    res.sheet = sheet_nr
        sheet_nr += 1

from dsspy.core import HelixType, HelixPositionType

def calculate_helices(residues):
    """
    Calculates helices from the H-bond information.
    """
    # Find helix patterns based on H-bonds
    for helix_type in [HelixType.ALPHA, HelixType._3_10, HelixType.PI]:
        stride = 3 + helix_type.value
        for i in range(len(residues) - stride):
            r1 = residues[i]
            r2 = residues[i + stride]
            if test_bond(r2, r1):
                if r1.helix_flags[helix_type] == HelixPositionType.END:
                    r1.helix_flags[helix_type] = HelixPositionType.START_AND_END
                else:
                    r1.helix_flags[helix_type] = HelixPositionType.START

                for j in range(1, stride):
                    res = residues[i+j]
                    if res.helix_flags[helix_type] == HelixPositionType.NONE:
                        res.helix_flags[helix_type] = HelixPositionType.MIDDLE

                r2.helix_flags[helix_type] = HelixPositionType.END

    # Assign secondary structure based on helix flags
    for i, res in enumerate(residues):
        # Alpha helices
        if res.helix_flags[HelixType.ALPHA] in [HelixPositionType.START, HelixPositionType.MIDDLE] and \
           i > 0 and residues[i-1].helix_flags[HelixType.ALPHA] in [HelixPositionType.START, HelixPositionType.MIDDLE]:
            for j in range(4):
                if i + j < len(residues):
                    residues[i+j].secondary_structure = StructureType.ALPHA_HELIX

        # 3-10 helices
        if res.helix_flags[HelixType._3_10] in [HelixPositionType.START, HelixPositionType.MIDDLE] and \
           i > 0 and residues[i-1].helix_flags[HelixType._3_10] in [HelixPositionType.START, HelixPositionType.MIDDLE]:
            is_empty = all(r.secondary_structure in [StructureType.LOOP, StructureType.HELIX_3] for r in residues[i:i+3])
            if is_empty:
                for j in range(3):
                    if i + j < len(residues):
                        residues[i+j].secondary_structure = StructureType.HELIX_3

        # Pi helices
        if res.helix_flags[HelixType.PI] in [HelixPositionType.START, HelixPositionType.MIDDLE] and \
           i > 0 and residues[i-1].helix_flags[HelixType.PI] in [HelixPositionType.START, HelixPositionType.MIDDLE]:
            is_empty = all(r.secondary_structure in [StructureType.LOOP, StructureType.HELIX_5] for r in residues[i:i+5])
            if is_empty:
                for j in range(5):
                    if i + j < len(residues):
                        residues[i+j].secondary_structure = StructureType.HELIX_5

    # Turns and Bends
    for i, res in enumerate(residues):
        if res.secondary_structure == StructureType.LOOP:
            is_turn = False
            for helix_type in [HelixType.ALPHA, HelixType._3_10, HelixType.PI]:
                stride = 3 + helix_type.value
                for k in range(1, stride):
                    if i >= k and residues[i-k].helix_flags[helix_type] == HelixPositionType.START:
                        is_turn = True
                        break
                if is_turn:
                    break

            if is_turn:
                res.secondary_structure = StructureType.TURN
            elif res.kappa > 70: # Bend condition from C++ code
                res.secondary_structure = StructureType.BEND


class DSSP:
    def __init__(self, residues):
        self.residues = residues

    def run(self):
        calculate_geometry(self.residues)
        calculate_hbonds(self.residues)
        calculate_beta_sheets(self.residues)
        calculate_helices(self.residues)
        # TODO: Add surface accessibility calculation here if time permits.

def dihedral_angle(p1, p2, p3, p4):
    """
    Calculates the dihedral angle between four points.
    Using the formula from https://en.wikipedia.org/wiki/Dihedral_angle
    """
    v12 = p1 - p2
    v43 = p4 - p3
    z = p2 - p3

    p = np.cross(z, v12)
    x = np.cross(z, v43)
    y = np.cross(z, x)

    u = np.dot(x, x)
    v = np.dot(y, y)

    if u > 0 and v > 0:
        u = np.dot(p, x) / np.sqrt(u)
        v = np.dot(p, y) / np.sqrt(v)
        if u != 0 or v != 0:
            return np.rad2deg(np.arctan2(v, u))
    return 360.0

def calculate_geometry(residues):
    """
    Calculates geometric properties of the residues, such as torsion angles.
    """
    for i, res in enumerate(residues):
        if res.prev_residue:
            # Phi angle (C-N-CA-C)
            try:
                p1 = res.prev_residue.biopython_residue['C'].get_coord()
                p2 = res.biopython_residue['N'].get_coord()
                p3 = res.biopython_residue['CA'].get_coord()
                p4 = res.biopython_residue['C'].get_coord()
                res.phi = dihedral_angle(p1, p2, p3, p4)
            except KeyError:
                pass # Atom not found

        if res.next_residue:
            # Psi angle (N-CA-C-N)
            try:
                p1 = res.biopython_residue['N'].get_coord()
                p2 = res.biopython_residue['CA'].get_coord()
                p3 = res.biopython_residue['C'].get_coord()
                p4 = res.next_residue.biopython_residue['N'].get_coord()
                res.psi = dihedral_angle(p1, p2, p3, p4)
            except KeyError:
                pass # Atom not found

            # Omega angle (CA-C-N-CA)
            try:
                p1 = res.biopython_residue['CA'].get_coord()
                p2 = res.biopython_residue['C'].get_coord()
                p3 = res.next_residue.biopython_residue['N'].get_coord()
                p4 = res.next_residue.biopython_residue['CA'].get_coord()
                res.omega = dihedral_angle(p1, p2, p3, p4)
            except KeyError:
                pass # Atom not found
