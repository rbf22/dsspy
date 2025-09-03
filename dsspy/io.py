"""
This module handles reading of PDB and mmCIF files.
"""

from Bio.PDB import PDBParser, MMCIFParser
from .core import Residue, ChainBreakType

def extract_residues(structure):
    """Extracts and numbers residues from a Biopython structure.

    This function iterates through a Biopython Structure object, creates a
    dsspy.core.Residue object for each valid amino acid residue, and sets
    up chain break information and links between adjacent residues.

    Args:
        structure (Bio.PDB.Structure.Structure): The Biopython structure object.

    Returns:
        tuple[list[dsspy.core.Residue], Bio.PDB.Structure.Structure]: A tuple
        containing the list of extracted Residue objects and the original
        structure object.
    """
    residues = []
    res_number = 0
    prev_res = None

    for model in structure:
        for chain in model:
            chain_break = ChainBreakType.NEW_CHAIN
            for residue in chain:
                # Skip heteroatoms
                if residue.get_id()[0] != ' ':
                    continue

                # In C++, there's a check for completeness, let's ensure the
                # backbone atoms are present
                if not ('N' in residue and 'CA' in residue and 'C' in residue and 'O' in residue):
                    continue

                res_number += 1
                dssp_res = Residue(residue, res_number)
                dssp_res.chain_break = chain_break
                chain_break = ChainBreakType.NONE # Reset after the first residue in a chain

                if prev_res:
                    # Check for chain gap
                    c_atom = prev_res.biopython_residue['C']
                    n_atom = dssp_res.biopython_residue['N']
                    # Peptide bond length is ~1.33A, 2.5A is a safe threshold
                    if c_atom - n_atom > 2.5:
                        dssp_res.chain_break = ChainBreakType.GAP
                        # In C++ code, the res number is incremented on a gap
                        res_number += 1

                    prev_res.next_residue = dssp_res
                    dssp_res.prev_residue = prev_res

                residues.append(dssp_res)
                prev_res = dssp_res

    return residues, structure

def read_pdb(pdb_file):
    """Reads a PDB file and extracts dsspy Residue objects.

    Args:
        pdb_file (file-like object): An open file handle to a PDB file.

    Returns:
        tuple[list[dsspy.core.Residue], Bio.PDB.Structure.Structure]: A tuple
        containing the list of extracted Residue objects and the parsed
        Biopython structure object.
    """
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    return extract_residues(structure)

def read_cif(cif_file):
    """Reads an mmCIF file and extracts dsspy Residue objects.

    Args:
        cif_file (file-like object): An open file handle to an mmCIF file.

    Returns:
        tuple[list[dsspy.core.Residue], Bio.PDB.Structure.Structure]: A tuple
        containing the list of extracted Residue objects and the parsed
        Biopython structure object.
    """
    parser = MMCIFParser()
    structure = parser.get_structure("protein", cif_file)
    return extract_residues(structure)
