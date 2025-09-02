import gzip
from io import StringIO
from dsspy.io import read_pdb, read_cif
from dsspy.core import Residue

def test_read_pdb():
    """
    Tests reading a PDB file.
    """
    # The original project has test files in the test/ directory.
    pdb_file_path = "test/reference_data/pdb1cbs.ent.gz"

    with gzip.open(pdb_file_path, "rt") as f:
        pdb_content = f.read()

    # The PDBParser can take a file-like object.
    pdb_file_like = StringIO(pdb_content)
    residues, _ = read_pdb(pdb_file_like)

    assert isinstance(residues, list)
    assert len(residues) > 0
    assert all(isinstance(r, Residue) for r in residues)


def test_read_cif():
    """
    Tests reading an mmCIF file.
    """
    cif_file_path = "test/reference_data/1cbs.cif.gz"

    with gzip.open(cif_file_path, "rt") as f:
        cif_content = f.read()

    cif_file_like = StringIO(cif_content)
    residues, _ = read_cif(cif_file_like)

    assert isinstance(residues, list)
    assert len(residues) > 0
    assert all(isinstance(r, Residue) for r in residues)
