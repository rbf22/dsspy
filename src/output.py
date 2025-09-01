import datetime
from .core import Residue, StructureType, HelixType, HelixPositionType

# 3-to-1 letter code for amino acids
AA_CODES = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'UNK': 'X'
}

def format_dssp_line(residue):
    """
    Formats a single residue into a line for the classic DSSP format.
    """
    res_num = residue.number
    pdb_seq_num = residue.biopython_residue.get_id()[1]
    pdb_ins_code = residue.biopython_residue.get_id()[2].strip() or ' '
    pdb_strand_id = residue.biopython_residue.get_full_id()[2]

    aa = AA_CODES.get(residue.resname, 'X')

    ss = residue.secondary_structure.value

    helix_flags = ''.join([
        '>' if residue.helix_flags[ht] == HelixPositionType.START else
        '<' if residue.helix_flags[ht] == HelixPositionType.END else
        'X' if residue.helix_flags[ht] == HelixPositionType.START_AND_END else
        str(ht.value + 3) if residue.helix_flags[ht] == HelixPositionType.MIDDLE and ht != HelixType.PP else
        'P' if residue.helix_flags[ht] == HelixPositionType.MIDDLE and ht == HelixType.PP else
        ' '
        for ht in [HelixType._3_10, HelixType.ALPHA, HelixType.PI, HelixType.PP]
    ])

    bend = 'S' if residue.bend else ' '

    chirality = '+' if residue.alpha > 0 and residue.alpha != 360 else '-' if residue.alpha < 0 else ' '

    bp1 = residue.beta_partner[0].number if residue.beta_partner[0] else 0
    bp2 = residue.beta_partner[1].number if residue.beta_partner[1] else 0

    # TODO: Implement bridge label and sheet label
    bridgelabel = ' '
    sheet = ' '

    acc = int(residue.accessibility)

    nho1 = f"{residue.hbond_acceptor[0]['residue'].number - res_num},{residue.hbond_acceptor[0]['energy']:.1f}" if residue.hbond_acceptor else "0, 0.0"
    onh1 = f"{residue.hbond_donor[0]['residue'].number - res_num},{residue.hbond_donor[0]['energy']:.1f}" if residue.hbond_donor else "0, 0.0"
    nho2 = f"{residue.hbond_acceptor[1]['residue'].number - res_num},{residue.hbond_acceptor[1]['energy']:.1f}" if len(residue.hbond_acceptor) > 1 else "0, 0.0"
    onh2 = f"{residue.hbond_donor[1]['residue'].number - res_num},{residue.hbond_donor[1]['energy']:.1f}" if len(residue.hbond_donor) > 1 else "0, 0.0"

    tco = f"{residue.tco:.3f}"
    kappa = f"{residue.kappa:.1f}"
    alpha = f"{residue.alpha:.1f}"
    phi = f"{residue.phi:.1f}"
    psi = f"{residue.psi:.1f}"

    x, y, z = residue.biopython_residue['CA'].get_coord()

    line = f"{res_num:>5d}{pdb_seq_num:>5d}{pdb_ins_code:>1}{pdb_strand_id:>1} {aa:>1}  {ss:>1}{helix_flags:>4}{bend:>1}{chirality:>1} {bp1:>4d}{bp2:>4d}{bridgelabel:>1}{sheet:>1} {acc:>4d} " \
           f"{nho1:>11s}{onh1:>11s}{nho2:>11s}{onh2:>11s}  " \
           f"{tco:>7s}{kappa:>6s}{alpha:>6s}{phi:>6s}{psi:>6s} " \
           f"{x:>6.1f}{y:>6.1f}{z:>6.1f}"

    return line

def write_dssp(residues, output_file):
    """
    Writes the DSSP output in the classic format to a file.
    """
    now = datetime.datetime.now()
    header = f"==== Secondary Structure Definition by the program DSSP, Python version ==== DATE={now.strftime('%Y-%m-%d')}        .\n" \
             "REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637\n"
    output_file.write(header)

    # TODO: Add more header lines from the C++ code (COMPND, SOURCE, AUTHOR, etc.)

    output_file.write("  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA\n")

    for residue in residues:
        output_file.write(format_dssp_line(residue) + '\n')
