"""
Tests for the io module.
"""
import gzip
from io import StringIO
from dsspy.io import read_pdb, read_cif
from dsspy.core import Residue, ChainBreakType

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

def test_read_pdb_incomplete_residue():
    """
    Tests that residues with missing backbone atoms are skipped.
    """
    pdb_content = """
ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   1.000   1.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.000   2.000   1.000  1.00  0.00           C
"""
    pdb_file_like = StringIO(pdb_content)
    residues, _ = read_pdb(pdb_file_like)
    assert len(residues) == 0

def test_read_pdb_chain_gap():
    """
    Tests that a gap in the chain is correctly identified.
    """
    pdb_content = """
ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   1.000   1.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.000   2.000   1.000  1.00  0.00           C
ATOM      4  O   ALA A   1       3.000   2.000   1.000  1.00  0.00           O
ATOM      5  N   ALA A   2      10.000  10.000  10.000  1.00  0.00           N
ATOM      6  CA  ALA A   2      11.000  10.000  10.000  1.00  0.00           C
ATOM      7  C   ALA A   2      11.000  11.000  10.000  1.00  0.00           C
ATOM      8  O   ALA A   2      12.000  11.000  10.000  1.00  0.00           O
"""
    pdb_file_like = StringIO(pdb_content)
    residues, _ = read_pdb(pdb_file_like)
    assert len(residues) == 2
    assert residues[1].chain_break == ChainBreakType.GAP
