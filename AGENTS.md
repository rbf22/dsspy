# Project Goals

The primary goal of this project is to convert the `mkdssp` C++ application into a modern Python package named `dsspy`. This conversion involves rewriting the C++ functionality in Python and leveraging standard libraries like `biopython` for handling protein structure data.

# Agent's Tasks

As the software engineering agent assigned to this project, my main responsibilities are:

1.  **Project Setup**: Establish a proper Python project structure using `pyproject.toml` for dependency management with `poetry`, and a `Makefile` for automating common development tasks.
2.  **Incremental Conversion**: Systematically convert the C++ code from the `libdssp` and `src` directories into Python code within the `dsspy` package.
3.  **Testing**: Write comprehensive unit tests for the new Python code to ensure its correctness and to maintain high code quality.
4.  **Dependency Management**: Utilize `biopython` for parsing PDB and mmCIF files, and manage all Python dependencies through `poetry`.

# Core Dependencies

The cornerstone of the Python conversion is the `biopython` library, which is used for all protein structure file parsing. All project dependencies are defined in the `pyproject.toml` file and managed by `poetry`.

# Makefile Usage

A `Makefile` is provided to simplify the development workflow. The key targets are:

-   `build`: Compiles the original C++ application using `cmake`. This is useful for reference and for testing the original behavior.
-   `poetry-install`: Installs the Python dependencies for the `dsspy` project into a virtual environment using `poetry`.
-   `run-tests`: Executes the Python unit test suite using `pytest`.
-   `run-linter`, `run-linter-fix`, `run-pylint`, `run-mypy`, `run-deptry`: These targets run a variety of linting and static analysis tools to ensure code quality.
-   `clean`: Removes the C++ build directory and other temporary artifacts.
