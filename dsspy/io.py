from Bio.PDB import PDBParser, MMCIFParser
from dsspy.core import Residue, ChainBreakType

def extract_residues(structure):
    """
    Extracts residues from a Biopython structure object and returns a list of dsspy.core.Residue objects.
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

                # In C++, there's a check for completeness, let's ensure the backbone atoms are present
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
                    if c_atom - n_atom > 2.5: # Peptide bond length is ~1.33A, 2.5A is a safe threshold
                        dssp_res.chain_break = ChainBreakType.GAP
                        res_number +=1 # In C++ code, the res number is incremented on a gap

                    prev_res.next_residue = dssp_res
                    dssp_res.prev_residue = prev_res

                residues.append(dssp_res)
                prev_res = dssp_res

    return residues

def read_pdb(pdb_file):
    """
    Reads a PDB file and returns a list of dsspy.core.Residue objects.
    """
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    return extract_residues(structure)

def read_cif(cif_file):
    """
    Reads an mmCIF file and returns a list of dsspy.core.Residue objects.
    """
    parser = MMCIFParser()
    structure = parser.get_structure("protein", cif_file)
    return extract_residues(structure)
