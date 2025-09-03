"""Command-line interface for dsspy.

This module provides a command-line tool to run DSSP calculations on PDB or
mmCIF files.
"""
import click
import sys
from .io import read_pdb, read_cif
from .hbond import calculate_h_bonds
from .secondary_structure import calculate_beta_sheets, calculate_pp_helices
from .accessibility import calculate_accessibility
from .output import write_dssp

@click.command()
@click.argument('input_file', type=click.Path(exists=True, dir_okay=False, readable=True))
@click.option('--output-file', '-o', type=click.Path(dir_okay=False, writable=True),
              help="Path to the output DSSP file (defaults to standard output).")
def main(input_file, output_file):
    """
    Assigns secondary structure to proteins using dsspy.

    Reads a PDB or mmCIF file and outputs a classic DSSP file.
    """
    try:
        # Determine file type and read the structure
        if input_file.lower().endswith('.cif') or input_file.lower().endswith('.mmcif'):
            with open(input_file, 'r', encoding='utf-8') as f:
                residues, structure = read_cif(f)
        elif input_file.lower().endswith('.pdb') or input_file.lower().endswith('.ent'):
            with open(input_file, 'r', encoding='utf-8') as f:
                residues, structure = read_pdb(f)
        else:
            click.echo("Error: Input file must be a PDB (.pdb, .ent) or mmCIF (.cif, .mmcif) file.", err=True)
            sys.exit(1)
    except Exception as e:
        click.echo(f"Error processing input file: {e}", err=True)
        sys.exit(1)

    # Run the DSSP pipeline
    calculate_h_bonds(residues)
    calculate_beta_sheets(residues)
    calculate_pp_helices(residues)
    calculate_accessibility(residues)

    # Write the output
    output_target = open(output_file, 'w', encoding='utf-8') if output_file else sys.stdout

    try:
        write_dssp(structure, residues, output_target)
    finally:
        if output_file:
            output_target.close()

if __name__ == '__main__':
    main()
