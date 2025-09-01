[![github CI](https://github.com/pdb-redo/dssp/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/pdb-redo/dssp/actions)
[![GitHub License](https://img.shields.io/github/license/pdb-redo/dssp)](https://github.com/pdb-redo/dssp/LICENSE)

DSSP 4.5
========

This is a rewrite of DSSP, now offering full mmCIF support. The difference
with previous releases of DSSP is that it now writes out an annotated mmCIF
file by default, storing the secondary structure information in the
`_struct_conf` category.

Another new feature in this version of DSSP is that it now defines
Poly-Proline helices as well.

The DSSP program was designed by _Wolfgang Kabsch_ and _Chris Sander_ to
standardize secondary structure assignment. DSSP is a database of secondary
structure assignments (and much more) for all protein entries in the Protein
Data Bank (PDB). DSSP is also the program that calculates DSSP entries from
PDB entries.

DSSP does **not** predict secondary structure.

Requirements
------------

A good, modern compiler is needed to build the mkdssp program since it uses
many new C++20 features.

Building
--------

To build the C++ version of `mkdssp`, you can use the `Makefile` in the root of the repository:

```console
make build
```

This will compile the code and place the artifacts in the `cpp_legacy/build` directory.

Alternatively, you can run the `cmake` commands directly from the root of the repository:

```console
cmake -S cpp_legacy -B cpp_legacy/build
cmake --build cpp_legacy/build
```

Python module
-------------

To build the Python module, you can use the `build` target in the `Makefile` at the root of the repository, which is configured to build the module by default.

```console
make build
```

Alternatively, you can run the `cmake` commands directly from the root of the repository:

```console
cmake -S cpp_legacy -B cpp_legacy/build -DBUILD_PYTHON_MODULE=ON
cmake --build cpp_legacy/build
```

After that you can use dssp in a python script, like this:

```python
from mkdssp import dssp
import os
import gzip

file_path = os.path.join("..", "test", "1cbs.cif.gz")

with gzip.open(file_path, "rt") as f:
    file_content = f.read()

dssp = dssp(file_content)

print("residues: ", dssp.statistics.residues)

for res in dssp:
    print(res.asym_id, res.seq_id, res.compound_id, res.type)

```

Usage
-----

See [manual page](doc/mkdssp.md) for more info. Or even better, see the [DSSP website](https://pdb-redo.eu/dssp).
