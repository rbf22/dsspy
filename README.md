# dsspy: A Python implementation of DSSP

![Build](https://github.com/rbf22/dsspy/actions/workflows/python-ci.yml/badge.svg)
![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)

`dsspy` is a Python package that provides a pure Python implementation of the DSSP algorithm for assigning secondary structure to proteins. This project is a complete rewrite of the original C++ `mkdssp` application, designed for seamless integration into modern Python-based bioinformatics workflows.

## Installation

This project uses [Poetry](https://python-poetry.org/) for dependency management and packaging. To install the necessary dependencies, follow these steps:

1.  **Install Poetry**:
    If you don't have Poetry installed, you can find instructions on the [official website](https://python-poetry.org/docs/#installation).

2.  **Install Dependencies**:
    Clone this repository and run the following command in the project root to install the required packages into a virtual environment:
    ```console
    poetry install
    ```

## Command-Line Usage

This package provides a command-line tool, `dsspy`, for running DSSP calculations directly from your terminal.

### Basic Usage

To run the calculation on a PDB or mmCIF file, simply provide the input file path. The DSSP output will be printed to standard output.

```console
poetry run dsspy test/reference_data/1cbs.cif
```

### Writing to a File

You can also specify an output file using the `-o` or `--output-file` option:

```console
poetry run dsspy test/reference_data/1cbs.cif -o 1cbs.dssp
```

### Getting Help

For a full list of options, use the `--help` flag:

```console
poetry run dsspy --help
```

## Library Usage

`dsspy` can also be used as a Python library. It processes protein structures in a sequential pipeline. Here is a complete example of how to read a structure from a CIF file, calculate all secondary structure features, and write the results to a classic DSSP file.


```python
import io
from dsspy.io import read_cif
from dsspy.hbond import calculate_h_bonds
from dsspy.secondary_structure import calculate_beta_sheets, calculate_pp_helices
from dsspy.accessibility import calculate_accessibility
from dsspy.output import write_dssp

# 1. Read a protein structure from a CIF file.
#    For this example, we use one of the test files.
with open("test/reference_data/1cbs.cif", "r") as f:
    residues, structure = read_cif(f)

# 2. Calculate H-bonds, which are necessary for the next steps.
calculate_h_bonds(residues)

# 3. Calculate secondary structures (beta sheets and helices).
calculate_beta_sheets(residues)
calculate_pp_helices(residues)

# 4. Calculate residue accessibility.
calculate_accessibility(residues)

# 5. Inspect the results.
#    The results of the calculations are stored as attributes on the Residue objects.
print("Residue | SS | Sheet | Accessibility")
print("--------|----|-------|---------------")
for res in residues[:5]:
    print(f"{res.number:>7} |  {res.secondary_structure.value} | {res.sheet:>5} | {res.accessibility:>10.1f}")

# 6. Write the output to a file in the classic DSSP format.
with open("output.dssp", "w") as f:
    write_dssp(structure, residues, f)

print("\nDSSP output written to output.dssp")
```

## Output Format

The output file (`output.dssp`) is in the classic DSSP format. Here is a snippet of what the output looks like:

```
  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
    1    1 A Cys  E              0   0  114    -2,-2.5    -1,-0.2    -2,-0.2     0, 0.0 -0.893 360.0-100.2-103.6 138.8   27.3   13.5   20.4
    2    2 A Ser  E              0   0   84    -2,-0.2     2,-2.5     0, 0.0     0, 0.0 -0.817 360.0-104.9 -93.8 143.8   30.5   14.2   21.9
    3    3 A Thr  E              0   0  113    -2,-2.5    -1,-0.2    -2,-0.2     0, 0.0 -0.933 360.0-104.9 -93.8 143.8   30.5   14.2   21.9
```

## Running Tests

To ensure the reliability and correctness of the `dsspy` package, a comprehensive test suite is provided. To run the tests, use the following command:

```console
poetry run pytest
```

This will execute the test suite and provide a coverage report.
