#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Protein Sequence Analysis Tool

This script analyzes protein sequences from a FASTA file to compute various
physicochemical properties. It utilizes the Biopython library for sequence
parsing and analysis, and pandas for data organization and output.

The calculated properties include:
- Sequence Length
- Molecular Weight (in kDa)
- Isoelectric Point (pI)
- Instability Index
- Grand Average of Hydropathicity (GRAVY)
- Aliphatic Index
- Net charge at pH 7

Results are compiled and can be saved into a CSV (default), TSV, or Excel (.xlsx) file.

The script is designed to be run from the command line, accepting input and
output file paths as arguments.

Dependencies:
- pandas
- biopython
- openpyxl (required by pandas for .xlsx export)

Author: Beaven Manjengwa
Date: September 30, 2025
Version: 1.0
"""

# Standard library imports
import argparse
import os
import sys

# Third-party library imports
try:
    import pandas as pd
    from Bio import SeqIO
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError as e:
    print(f"Error: A required library is missing. {e}", file=sys.stderr)
    print("Please install the necessary packages using: pip install -r requirements.txt", file=sys.stderr)
    sys.exit(1)


def analyze_protein_sequences(fasta_file: str, output_file: str) -> None:
    """
    Analyzes protein sequences from a FASTA file, calculates key
    physicochemical properties, and saves the results to a file.

    The output format is determined by the extension of the output_file name
    (e.g., .csv, .tsv, .xlsx).

    Args:
        fasta_file (str): The full path to the input FASTA file containing
                          protein sequences.
        output_file (str): The full path where the output file will be saved.
    
    Returns:
        None
    """
    # List to hold the dictionaries of results for each protein sequence.
    data = []

    try:
        # Open and parse the FASTA file.
        with open(fasta_file, 'r') as file:
            for record in SeqIO.parse(file, "fasta"):
                # Ensure the sequence is a string for analysis.
                # Sanitize sequence by removing stop codons ('*') and whitespace.
                sequence = str(record.seq).upper().replace('*', '').strip()
                
                # Skip empty sequences which would cause errors
                if not sequence:
                    print(f"Warning: Skipping empty sequence for ID: {record.id}")
                    continue

                # Create a ProteinAnalysis object from the sequence.
                analysed_seq = ProteinAnalysis(sequence)

                # --- Property Calculations ---
                
                # Get amino acid composition using the modern attribute to avoid deprecation warning.
                aa_percent = analysed_seq.amino_acids_percent

                # Manual calculation of the Aliphatic Index for backward compatibility.
                # Formula (Ikai, 1980): Aliphatic index = X(Ala) + a * X(Val) + b * (X(Ile) + X(Leu))
                # Note: amino_acids_percent provides percentages directly, so no need to multiply by 100.
                aliphatic_index = (aa_percent.get('A', 0)) + \
                                  (2.9 * aa_percent.get('V', 0)) + \
                                  (3.9 * (aa_percent.get('I', 0) + aa_percent.get('L', 0)))

                # Compile all results into a dictionary for this sequence.
                result = {
                    "gene_id": record.id,
                    "sequence_length": len(sequence),
                    "molecular_weight_kDa": round(analysed_seq.molecular_weight() / 1000, 2),
                    "isoelectric_point_pI": round(analysed_seq.isoelectric_point(), 2),
                    "instability_index": round(analysed_seq.instability_index(), 2),
                    "gravy": round(analysed_seq.gravy(), 3),
                    "aliphatic_index": round(aliphatic_index, 2),
                    "charge_at_pH7": round(analysed_seq.charge_at_pH(7.0), 2)
                }

                # Add the result dictionary to our main data list.
                data.append(result)

        # After processing all sequences, convert the list of dictionaries to a pandas DataFrame.
        if not data:
            print("Warning: No valid sequences were found in the input file. Output file will not be created.")
            return

        df = pd.DataFrame(data)
        
        # Determine the output format from the file extension.
        output_format = os.path.splitext(output_file)[1].lower()

        if output_format == '.csv':
            df.to_csv(output_file, index=False)
        elif output_format == '.tsv':
            df.to_csv(output_file, index=False, sep='\t')
        elif output_format == '.xlsx':
            df.to_excel(output_file, index=False, engine='openpyxl')
        else:
            print(f"Error: Unsupported output file format '{output_format}'. Please use .csv, .tsv, or .xlsx.", file=sys.stderr)
            sys.exit(1)

        # Provide user feedback upon successful completion.
        print(f"\nAnalysis complete! Processed {len(data)} sequences.")
        print(f"Results have been saved to: {output_file}")

    except FileNotFoundError:
        print(f"Error: The input file '{fasta_file}' was not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        # Catch any other potential errors during file processing or analysis.
        print(f"An unexpected error occurred during analysis: {str(e)}", file=sys.stderr)
        sys.exit(1)


def main() -> None:
    """
    Main function to handle command-line argument parsing and script execution flow.
    """
    # Initialize the argument parser with a description and usage examples.
    parser = argparse.ArgumentParser(
        description="Analyze protein sequences from a FASTA file and calculate various properties.",
        formatter_class=argparse.RawTextHelpFormatter, # Using RawTextHelpFormatter for better control over epilog formatting
        epilog="""
Examples of usage:
  # Basic usage with a default output file (e.g., input.properties.csv)
  python protein_analysis.py -i input.fasta

  # Specify a CSV output file
  python protein_analysis.py -i input.fasta -o results.csv

  # Specify a TSV output file
  python protein_analysis.py --input sequences.fasta --output results.tsv

  # Specify an Excel output file and create a new directory for it
  python protein_analysis.py --input data/proteins.fasta --output results/analysis.xlsx
        """
    )

    # --- Argument Definitions ---
    # CHANGE: Switched to required, flagged arguments for clarity. This is more standard for CLI tools.
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Path to the input FASTA file.'
    )
    parser.add_argument(
        '-o', '--output',
        help='Path for the output file. Format (.csv, .tsv, .xlsx) is determined by extension. '
             'If omitted, a .csv file is created based on the input filename.'
    )

    # Parse the arguments provided by the user.
    args = parser.parse_args()

    # --- Argument Validation and Processing ---
    input_file = args.input
    output_file = args.output

    # If no output file is specified, auto-generate a .csv filename based on the input.
    if not output_file:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}.properties.csv"
        print(f"No output file specified. Defaulting to: {output_file}")

    # Ensure the directory for the output file exists. If not, create it.
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
        except OSError as e:
            print(f"Error: Could not create output directory '{output_dir}'. {e}", file=sys.stderr)
            sys.exit(1)

    # Call the main analysis function with the validated file paths.
    analyze_protein_sequences(input_file, output_file)


# Standard Python entry point.
# Ensures that the main() function is called only when the script is executed directly.
if __name__ == "__main__":
    main()
