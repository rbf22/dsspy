# Future Refactoring Tasks

This file lists the major refactoring tasks that still need to be completed to fully convert the C++ `mkdssp` application to the Python `dsspy` package.

### 1. Solvent Accessibility Calculation

- **Task:** Implement the solvent accessibility calculation in `dsspy/algorithm.py`.
- **Details:** The current Python implementation is missing this feature entirely. The C++ code (`libdssp/src/dssp.cpp`) contains a function `residue::CalculateSurface` that implements this. The algorithm involves creating a sphere of points around each atom (using a Fibonacci sphere) and checking how many of those points are exposed to the solvent.
- **Action:** Port the `CalculateSurface` and `CalculateAccessibilities` logic to Python. This will likely involve creating a new function `calculate_accessibility` in `dsspy/algorithm.py` and updating the `Residue` class in `dsspy/core.py` to store the accessibility value.

### 2. DSSP Header Generation

- **Task:** Implement the generation of the full DSSP header in `dsspy/output.py`.
- **Details:** The current output has a minimal header. The C++ implementation extracts detailed information from the input file (e.g., `COMPND`, `SOURCE`, `AUTHOR` records).
- **Action:** Add functions to `dsspy/output.py` that extract the relevant metadata from the `biopython` structure object and format it into the classic DSSP header format.

### 3. Detailed Output Formatting

- **Task:** Complete the formatting of all fields in the DSSP output line in `dsspy/output.py`.
- **Details:** The current output formatting has `TODO`s for the `bridgelabel` and `sheet` columns.
- **Action:**
    - The `sheet` number is already calculated; it just needs to be correctly formatted in the output string.
    - The `bridgelabel` requires storing and retrieving bridge partner information in the `Residue` object. This information needs to be passed to the output formatting function.

### 4. Poly-proline (PPII) Helix Calculation

- **Task:** Review and verify the poly-proline helix calculation.
- **Details:** The C++ code has a function `CalculatePPHelices`. The Python implementation has a placeholder for `HelixType.PP` but does not seem to have the full logic for identifying these helices based on phi/psi angles.
- **Action:** Port the logic from `CalculatePPHelices` in `dssp.cpp` to `dsspy/algorithm.py` and update the `Residue` and `HelixType` data structures as needed.
