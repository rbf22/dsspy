import re
from click.testing import CliRunner
from dsspy.cli import main
import os

def test_cli_with_cif_input_to_stdout():
    """
    Tests that the CLI runs with a CIF input and prints to stdout.
    """
    runner = CliRunner()
    input_file = "test/reference_data/1cbs.cif"

    result = runner.invoke(main, [input_file])

    assert result.exit_code == 0
    assert result.output
    assert "==== Secondary Structure Definition by the program DSSP" in result.output
    # Check for a line that starts with residue 5 and has structure E
    assert re.search(r"^\s*5\s+5\s+A\s+G\s+E", result.output, re.MULTILINE)

def test_cli_to_file():
    """
    Tests that the CLI can write to an output file.
    """
    runner = CliRunner()
    input_file_path = "test/reference_data/1cbs.cif"
    output_file_name = "test_output.dssp"

    with open(input_file_path, 'r', encoding='utf-8') as f:
        cif_content = f.read()

    # Use the runner's isolated filesystem for the output file
    with runner.isolated_filesystem():
        # Create a dummy input file in the isolated filesystem
        with open("input.cif", "w", encoding='utf-8') as f:
            f.write(cif_content)

        result = runner.invoke(main, ["input.cif", '-o', output_file_name])

        assert result.exit_code == 0
        assert os.path.exists(output_file_name)
        with open(output_file_name, 'r', encoding='utf-8') as f:
            content = f.read()
            assert "==== Secondary Structure Definition by the program DSSP" in content
            # Check for a line that starts with residue 5 and has structure E
            assert re.search(r"^\s*5\s+5\s+A\s+G\s+E", content, re.MULTILINE)

def test_cli_file_not_found():
    """
    Tests that the CLI exits with an error if the input file is not found.
    """
    runner = CliRunner()
    input_file = "non_existent_file.pdb"

    result = runner.invoke(main, [input_file])

    assert result.exit_code != 0
    # Click provides its own error message for non-existent paths
    assert "Error: Invalid value for 'INPUT_FILE'" in result.output
