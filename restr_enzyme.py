#!/usr/bin/env python3
"""
Protein Digestion Tool

This script reads a Fasta formatted file, digests protein sequences using specified enzymes,
and writes the resulting peptides to an output file.
"""

import argparse
import os

# Fasta file reading function
def read_fasta(filename):
    """Reads and parses a Fasta formatted file."""
    try:
        with open(filename, "r") as fastafile:
            lines = fastafile.readlines()
    except FileNotFoundError:
        print("Error: File not found.")
        return [], {}

    order = []       # List to store sequence names in the order they appear
    sequences = {}   # Dictionary to store sequences keyed by sequence names
    current_sequence = ""  # Temporary variable to build the current sequence
    current_name = ""      # Temporary variable to store the name of the current sequence

    for line in lines:
        if line.startswith('>'):  # Header line indicating a new sequence
            if current_sequence:
                sequences[current_name] = current_sequence
            current_name = line.split()[0]  # Extract sequence name
            order.append(current_name)      # Add sequence name to the order list
            current_sequence = ""           # Reset current sequence
        else:  # Sequence line
            current_sequence += line.strip()  # Add sequence without leading/trailing whitespaces

    # Add the last sequence after the loop ends
    if current_sequence:
        sequences[current_name] = current_sequence

    return order, sequences

# Digestion function
def digest(sequence, enzyme):
    """Digests a protein sequence using specified enzyme."""
    peptides = []    # List to store generated peptides
    temp_peptide = ""  # Temporary variable to build the current peptide
    cut_list = []     # List of amino acids where the enzyme cuts

    # Define the cut list based on the selected enzyme
    if enzyme == 't':  # Trypsin
        cut_list = ['K', 'R']
    elif enzyme == 'r':  # Endoproteinase Arg-C
        cut_list = ['R']
    elif enzyme == 'k':  # Endoproteinase Lys-C
        cut_list = ['K']
    elif enzyme == 'e':  # V8 proteinase (Glu-C)
        cut_list = ['E']

    # Iterate over each amino acid in the sequence
    for aa in sequence:
        temp_peptide += aa  # Add the amino acid to the temporary peptide sequence
        # Check if the current amino acid is a cutting point
        if aa in cut_list:
            if len(temp_peptide) > 1:  # Ensure the peptide is not empty
                peptides.append(temp_peptide)  # Add the peptide to the list
            temp_peptide = ""  # Reset the temporary peptide sequence

    # Add the last peptide after the loop ends
    if temp_peptide:
        peptides.append(temp_peptide)

    return peptides

if __name__ == "__main__":
    # Dictionary mapping enzyme codes to enzyme names
    enzymes = {"t": "trypsin", "r": "endoproteinase Arg-C",
               "k": "endoproteinase Lys-C", "e": "V8 proteinase (Glu-C)"}
    # List of valid minimum peptide lengths
    pepmins = list(range(1, 11))

    # Set up argument parser
    parser = argparse.ArgumentParser(description='Digest proteins from a Fasta file')
    parser.add_argument("input_file", help="input Fasta file")
    parser.add_argument("-e", "--enzyme", choices=enzymes.keys(), default='t',
                        help="enzyme to use for digestion (default: trypsin)")
    parser.add_argument("-p", "--pepmin", type=int, choices=pepmins, default=5,
                        help="minimum peptide length to include in the results (default: 5)")
    args = parser.parse_args()

    # Extract arguments
    input_file = args.input_file
    enzyme = args.enzyme
    pepmin = args.pepmin

    # Check if input file exists
    if not os.path.isfile(input_file):
        print("Error: Input file not found.")
        exit(1)

    # Read sequences from input file
    order, sequences = read_fasta(input_file)

    # Check if any sequences were read
    if not order:
        print("Error: No sequences found in the input file.")
        exit(1)

    # Set up output file name
    output_file = f"output_{enzymes[enzyme]}.txt"

    # Digest sequences and write peptides to output file
    with open(output_file, 'w') as outfile:
        for prot in order:
            peptides = digest(sequences[prot], enzyme)
            for i, peptide in enumerate(peptides):
                if len(peptide) >= pepmin:
                    outfile.write(f"{prot} peptide {i+1:03d}\n{peptide}\n")
                    print(f"{prot} peptide {i+1:03d}\n{peptide}")

   # print("Digestion completed. Results written to", output_file)
