"""
This module handles the output of DSSP files in the classic format.
"""

import datetime
from .core import HelixType, HelixPositionType
from .utils import get_sheet_label

# 3-to-1 letter code for amino acids
AA_CODES = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'UNK': 'X'
}

def format_dssp_line(residue):
    """Formats a single residue into a line for the classic DSSP output format.

    This function takes a single Residue object, which must have all its
    DSSP parameters calculated, and formats it into a single line string
    that conforms to the classic DSSP file format.

    Args:
        residue (dsspy.core.Residue): The residue to format.

    Returns:
        str: A string representing the residue in classic DSSP format.
    """
    # pylint: disable=too-many-locals
    res_num = residue.number
    pdb_seq_num = residue.biopython_residue.get_id()[1]
    pdb_ins_code = residue.biopython_residue.get_id()[2].strip() or ' '
    pdb_strand_id = residue.biopython_residue.get_full_id()[2]

    aa = AA_CODES.get(residue.resname, 'X')

    ss = residue.secondary_structure.value

    helix_flags_list = []
    for ht in [HelixType.THREE_TEN, HelixType.ALPHA, HelixType.PI, HelixType.PP]:
        flag = residue.helix_flags[ht]
        if flag == HelixPositionType.START:
            helix_flags_list.append('>')
        elif flag == HelixPositionType.END:
            helix_flags_list.append('<')
        elif flag == HelixPositionType.START_AND_END:
            helix_flags_list.append('X')
        elif flag == HelixPositionType.MIDDLE:
            if ht == HelixType.PP:
                helix_flags_list.append('P')
            else:
                helix_flags_list.append(str(ht.value + 3))
        else:
            helix_flags_list.append(' ')
    helix_flags = "".join(helix_flags_list)

    bend = 'S' if residue.bend else ' '

    chirality = '+' if residue.alpha > 0 and residue.alpha != 360 else '-' if residue.alpha < 0 else ' '

    bp1 = residue.beta_partner[0].residue.number if residue.beta_partner[0].residue else 0
    bp2 = residue.beta_partner[1].residue.number if residue.beta_partner[1].residue else 0

    sheet = get_sheet_label(residue.sheet)

    bridgelabel = ' '
    if residue.beta_partner[0].residue:
        ladder = residue.beta_partner[0].ladder
        if ladder is not None:
            bridgelabel = chr(ord('a') + ladder % 26)

    acc = int(residue.accessibility)

    nho1 = (f"{residue.hbond_acceptor[0].residue.number - res_num},"
            f"{residue.hbond_acceptor[0].energy:.1f}"
            if residue.hbond_acceptor[0].residue else "0, 0.0")
    onh1 = (f"{residue.hbond_donor[0].residue.number - res_num},"
            f"{residue.hbond_donor[0].energy:.1f}"
            if residue.hbond_donor[0].residue else "0, 0.0")
    nho2 = (f"{residue.hbond_acceptor[1].residue.number - res_num},"
            f"{residue.hbond_acceptor[1].energy:.1f}"
            if residue.hbond_acceptor[1].residue else "0, 0.0")
    onh2 = (f"{residue.hbond_donor[1].residue.number - res_num},"
            f"{residue.hbond_donor[1].energy:.1f}"
            if residue.hbond_donor[1].residue else "0, 0.0")

    tco = f"{residue.tco:.3f}"
    kappa = f"{residue.kappa:.1f}"
    alpha = f"{residue.alpha:.1f}"
    phi = f"{residue.phi:.1f}"
    psi = f"{residue.psi:.1f}"

    x, y, z = residue.biopython_residue['CA'].get_coord()

    line = (  # pylint: disable=line-too-long
        f"{res_num:>5d}{pdb_seq_num:>5d}{pdb_ins_code:>1}{pdb_strand_id:>1} {aa:>1}  "
        f"{ss:>1}{helix_flags:>4}{bend:>1}{chirality:>1} {bp1:>4d}{bp2:>4d}"
        f"{bridgelabel:>1}{sheet:>1} {acc:>4d} "
        f"{nho1:>11s}{onh1:>11s}{nho2:>11s}{onh2:>11s}  "
        f"{tco:>7s}{kappa:>6s}{alpha:>6s}{phi:>6s}{psi:>6s} "
        f"{x:>6.1f}{y:>6.1f}{z:>6.1f}"
    )

    return line

def _format_header(header_dict):
    """
    Formats the header dictionary into the classic DSSP header format.
    """
    lines = []
    now = datetime.datetime.now()

    lines.append("==== Secondary Structure Definition by the program DSSP, "
                 f"Python version ==== DATE={now.strftime('%Y-%m-%d')}        .")
    lines.append("REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637")

    head = header_dict.get('head', '').upper()
    dep_date = header_dict.get('deposition_date', '')
    id_code = header_dict.get('idcode', '')
    header_line = f"HEADER    {head:<40s} {dep_date:<9s}   {id_code:<4s}"
    lines.append(header_line)

    compnd_info = []
    # The 'compound' entry in biopython is a list of dictionaries.
    for comp in header_dict.get('compound', []):
        for key, value in comp.items():
            compnd_info.append(f"{key.upper()}: {value}")
    lines.append("COMPND    " + "; ".join(compnd_info))

    source_info = []
    for src in header_dict.get('source', []):
        for key, value in src.items():
            source_info.append(f"{key.upper()}: {value}")
    lines.append("SOURCE    " + "; ".join(source_info))

    author_line = "AUTHOR    " + header_dict.get('author', '')
    lines.append(author_line)

    # TODO: Add statistics section

    return "\n".join(lines)


def write_dssp(structure, residues, output_file):
    """Writes the full DSSP output in the classic format to a file.

    This function generates the complete DSSP output, including the header
    and the per-residue data lines, and writes it to the specified file-like
    object.

    Args:
        structure (Bio.PDB.Structure.Structure): The Biopython structure object,
            used to generate the header.
        residues (list[dsspy.core.Residue]): The list of Residue objects with
            all DSSP parameters calculated.
        output_file (file-like object): An open, writable file handle to write
            the output to.
    """
    header_text = _format_header(structure.header)
    output_file.write(header_text + '\n')

    output_file.write("  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    "  # pylint: disable=line-too-long
                      "N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA\n")

    for residue in residues:
        output_file.write(format_dssp_line(residue) + '\n')
