# dsspy: A Python implementation of DSSP

![Build](https://github.com/rbf22/dsspy/actions/workflows/python-ci.yml/badge.svg)
![License](https://img.shields.io/github/license/rbf22/dsspy)

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

## Usage

Here's a basic example of how to use `dsspy` to analyze a protein structure from a PDB file:

```python
from dsspy.io import DSSP

# Example of reading a PDB file and calculating DSSP
with open("path/to/your/protein.pdb", "r") as f:
    dssps = DSSP.from_pdb(f)

# Accessing DSSP data
for d in dssps:
    print(f"Residue {d.resnum}: {d.aa} - {d.ss}")

```

## Running Tests

To ensure the reliability and correctness of the `dsspy` package, a comprehensive test suite is provided. To run the tests, use the following command:

```console
poetry run pytest
```

This will execute the test suite and provide a coverage report.

## Legacy C++ Version

The original C++ version of `mkdssp` can be found on the pdbredo github page.
