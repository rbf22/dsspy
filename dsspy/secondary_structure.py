"""
This module contains functions for calculating secondary structure elements
like beta sheets and helices.
"""
from .core import (
    Residue,
    Bridge,
    BridgeType,
    StructureType,
    BridgePartner,
    HelixPositionType,
    HelixType,
)


def calculate_beta_sheets(residues: list[Residue]):
    """Calculates and assigns beta sheet information based on H-bond patterns.

    This function identifies beta bridges, extends them into ladders, and groups
    ladders into sheets. The results are assigned back to the `Residue` objects
    in the input list. This is a core part of the DSSP algorithm.

    Note: This function has side effects, as it modifies the Residue objects
    in the input list.

    Args:
        residues (list[Residue]): A list of all Residue objects in the structure.
            These objects must have their H-bond information pre-calculated.
    """
    # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    max_hbond_energy = -0.5

    def _test_bond(res1, res2):
        """Checks if res1 has an H-bond to res2."""
        if res1 is None or res2 is None:
            return False
        for hbond in res1.hbond_acceptor:
            if hbond.residue == res2 and hbond.energy < max_hbond_energy:
                return True
        return False

    def _no_chain_break(res1, res2):
        """Checks if there is a continuous chain of residues between res1 and res2."""
        if res1 is None or res2 is None:
            return False
        if res1.biopython_residue.get_parent() != res2.biopython_residue.get_parent():
            return False
        if res1.number > res2.number:
            res1, res2 = res2, res1

        current = res1
        while current and current.number < res2.number:
            if current.next_residue is None or current.next_residue.number != current.number + 1:
                return False
            current = current.next_residue
        return current is not None and current.number == res2.number

    def _test_bridge(res1, res2):
        """
        Tests for a parallel or anti-parallel bridge between res1 and res2.
        This is a port of the `TestBridge` function from the C++ dssp implementation.
        """
        a = res1.prev_residue
        b = res1
        c = res1.next_residue
        d = res2.prev_residue
        e = res2
        f = res2.next_residue

        if not (a and c and _no_chain_break(a, c) and d and f and _no_chain_break(d, f)):
            return BridgeType.NONE

        # Pattern for parallel bridge
        if (_test_bond(a, e) and _test_bond(e, c)) or (_test_bond(d, b) and _test_bond(b, f)):
            return BridgeType.PARALLEL

        # Pattern for anti-parallel bridge
        if (_test_bond(b, e) and _test_bond(e, b)) or (_test_bond(a, f) and _test_bond(d, c)):
            return BridgeType.ANTIPARALLEL

        return BridgeType.NONE

    # 1. Find initial bridges
    bridges = []
    for i, res1 in enumerate(residues):
        for j in range(i + 1, len(residues)):
            res2 = residues[j]

            if res1.biopython_residue.get_parent() != res2.biopython_residue.get_parent():
                continue

            # Residues in a bridge must be separated by at least 2
            if abs(res1.number - res2.number) < 3:
                continue

            bridge_type = _test_bridge(res1, res2)
            if bridge_type != BridgeType.NONE:
                bridges.append(Bridge(res1, res2, bridge_type))

    # print(f"Found {len(bridges)} initial bridges")

    # 2. Extend ladders (with bulge logic)
    bridges.sort(key=lambda b: (b.i[0].biopython_residue.get_parent().id, b.i[0].number))

    i = 0
    while i < len(bridges):
        j = i + 1
        while j < len(bridges):
            b_i = bridges[i]
            b_j = bridges[j]

            merged = False

            if b_i.type == b_j.type and \
               _no_chain_break(b_i.i[0], b_j.i[-1]) and \
               _no_chain_break(b_i.j[0], b_j.j[-1]):

                iei = b_i.i[-1].number
                ibj = b_j.i[0].number

                if b_i.type == BridgeType.PARALLEL:
                    jei = b_i.j[-1].number
                    jbj = b_j.j[0].number
                    if (jbj - jei < 6 and ibj - iei < 3) or (jbj - jei < 3):
                        b_i.i.extend(b_j.i)
                        b_i.j.extend(b_j.j)
                        merged = True
                else:  # Antiparallel
                    jbi = b_i.j[0].number
                    jej = b_j.j[-1].number
                    if (jbi - jej < 6 and ibj - iei < 3) or (jbi - jej < 3):
                        b_i.i.extend(b_j.i)
                        # For antiparallel, the j-strand of the second bridge is added to the front
                        b_i.j = b_j.j + b_i.j
                        merged = True

            if merged:
                bridges.pop(j)
            else:
                j += 1
        i += 1

    # print(f"Found {len(bridges)} ladders after extension")

    # 3. Create sheets from linked ladders
    ladderset = set(bridges)
    sheet_nr = 1
    ladder_nr = 0
    while ladderset:
        sheet = {ladderset.pop()}

        while True:
            newly_added = set()
            for b1 in sheet:
                for b2 in ladderset:
                    # Check if any residue from b2 is a neighbor of any residue in b1
                    b1_res = set(b1.i) | set(b1.j)
                    b2_res = set(b2.i) | set(b2.j)
                    is_linked = any(
                        r1.next_residue in b2_res or r1.prev_residue in b2_res
                        for r1 in b1_res
                    )
                    if is_linked:
                        newly_added.add(b2)

            if not newly_added:
                break

            sheet.update(newly_added)
            ladderset.difference_update(newly_added)

        for bridge in sheet:
            bridge.sheet = sheet_nr
            bridge.ladder = ladder_nr
            ladder_nr += 1
        sheet_nr += 1

    # 4. Annotate residues with bridge and sheet info
    res_map = {res.number: res for res in residues}
    for bridge in bridges:
        ss = StructureType.BETA_BRIDGE if len(bridge.i) == 1 else StructureType.STRAND

        # Find which beta partner slot to fill (0 or 1)
        beta_i = 1 if bridge.i[0].beta_partner[0].residue is not None else 0
        beta_j = 1 if bridge.j[0].beta_partner[0].residue is not None else 0

        if bridge.type == BridgeType.PARALLEL:
            for k, res_i in enumerate(bridge.i):
                res_j = bridge.j[k]
                res_i.beta_partner[beta_i] = BridgePartner(res_j, bridge.ladder, True)
                res_j.beta_partner[beta_j] = BridgePartner(res_i, bridge.ladder, True)
        else:  # Antiparallel
            for k, res_i in enumerate(bridge.i):
                res_j = bridge.j[-(k + 1)]
                res_i.beta_partner[beta_i] = BridgePartner(res_j, bridge.ladder, False)
                res_j.beta_partner[beta_j] = BridgePartner(res_i, bridge.ladder, False)

        # Annotate all residues in the strand range
        i_start = bridge.i[0].number
        i_end = bridge.i[-1].number
        for res_num in range(i_start, i_end + 1):
            res = res_map.get(res_num)
            if res:
                if res.secondary_structure != StructureType.STRAND:
                    res.secondary_structure = ss
                res.sheet = bridge.sheet

        j_start = min(r.number for r in bridge.j)
        j_end = max(r.number for r in bridge.j)
        for res_num in range(j_start, j_end + 1):
            res = res_map.get(res_num)
            if res:
                if res.secondary_structure != StructureType.STRAND:
                    res.secondary_structure = ss
                res.sheet = bridge.sheet

def calculate_pp_helices(residues: list[Residue], stretch_length: int = 3):
    """Identifies Polyproline II (PPII) helices based on phi/psi angles.

    This function scans through the residues and identifies stretches that have
    the characteristic phi/psi angles of a Polyproline II helix.

    Note: This function has side effects, as it modifies the Residue objects
    in the input list.

    Args:
        residues (list[Residue]): A list of all Residue objects in the structure.
        stretch_length (int, optional): The minimum number of consecutive
            residues that must match the PPII criteria to be considered a helix.
            Defaults to 3.
    """
    n_residues = len(residues)
    if n_residues < stretch_length:
        return

    epsilon = 29.0
    phi_min = -75.0 - epsilon
    phi_max = -75.0 + epsilon
    psi_min = 145.0 - epsilon
    psi_max = 145.0 + epsilon

    for i in range(n_residues - stretch_length + 1):
        in_helix = True
        for j in range(stretch_length):
            res = residues[i + j]
            if not (res.phi and res.psi and
                    phi_min <= res.phi <= phi_max and
                    psi_min <= res.psi <= psi_max):
                in_helix = False
                break

        if in_helix:
            for j in range(stretch_length):
                res = residues[i + j]
                if res.secondary_structure == StructureType.LOOP:
                    res.secondary_structure = StructureType.HELIX_PPII

                if j == 0: # Start of the helix
                    if res.helix_flags[HelixType.PP] == HelixPositionType.NONE:
                        res.helix_flags[HelixType.PP] = HelixPositionType.START
                    elif res.helix_flags[HelixType.PP] == HelixPositionType.END:
                        res.helix_flags[HelixType.PP] = HelixPositionType.START_AND_END
                elif j == stretch_length - 1:  # End of the helix
                    res.helix_flags[HelixType.PP] = HelixPositionType.END
                else: # Middle of the helix
                    res.helix_flags[HelixType.PP] = HelixPositionType.MIDDLE
