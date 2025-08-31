from Bio.PDB import PDBParser, MMCIFParser

def read_pdb(pdb_file):
    """
    Reads a PDB file and returns a structure object.
    """
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)
    return structure

def read_cif(cif_file):
    """
    Reads an mmCIF file and returns a structure object.
    """
    parser = MMCIFParser()
    structure = parser.get_structure("protein", cif_file)
    return structure
