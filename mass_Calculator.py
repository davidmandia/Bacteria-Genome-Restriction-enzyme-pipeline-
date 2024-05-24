#!/usr/bin/env python3
"""
Peptide Mass Calculation Tool

This script reads peptide sequences from a Fasta file, calculates their masses,
and writes the results to an output file.
"""

import argparse

# Function to read Fasta formatted file
def fastaread(filename):
    """Reads and parses a Fasta formatted file."""
    order = []   # List to store sequence names in the order they appear
    seqs = {}    # Dictionary to store sequences keyed by sequence names

    with open(filename, "r") as fastafile:
        for line in fastafile:
            if line.startswith('>'):
                seq_name = line.split()[0][1:] + ',' + line.split()[2]  # Extract sequence name
                order.append(seq_name)
                seqs[seq_name] = ''
            else:
                seqs[seq_name] = line.strip()

    return order, seqs

# Function to calculate peptide mass
def pep2mass(seq, type):
    """Calculates the mass of a peptide sequence."""
    # Dictionary mapping amino acids to their monoisotopic or average masses
    mono = {'A': 71.0371, 'C': 103.0092, 'D': 115.0269, 'E': 129.0426, 'F': 147.0684,
            'G': 57.0215, 'H': 137.0589, 'I': 113.0841, 'K': 128.0950, 'L': 113.0841,
            'M': 131.0405, 'N': 114.0429, 'P': 97.0528, 'Q': 128.0586, 'R': 156.1011,
            'S': 87.0320, 'T': 101.0477, 'V': 99.0684, 'W': 186.0793, 'Y': 163.0633,
            '*': 0.0}  # Stop codon mass is 0

    aver = {'A': 71.08, 'C': 103.14, 'D': 115.09, 'E': 129.12, 'F': 147.18,
            'G': 57.05, 'H': 137.14, 'I': 113.16, 'K': 128.17, 'L': 113.16,
            'M': 131.19, 'N': 114.10, 'P': 97.12, 'Q': 128.13, 'R': 156.19,
            'S': 87.08, 'T': 101.10, 'V': 99.13, 'W': 186.21, 'Y': 163.18,
            '*': 0.0}

    # Define water mass based on the type of mass calculation
    water = 19.0153 if type == "a" else 19.0106

    # Select mass dictionary based on the type of mass calculation
    mass_dict = mono if type == "mono" else aver

    # Calculate mass of the peptide
    mass = round(sum(mass_dict[aa] for aa in seq) + water, 4)

    return mass

if __name__ == "__main__":
    # Define valid arguments for mass type
    masstypes = ['a', 'm']

    # Set up argument parser
    parser = argparse.ArgumentParser(description='Convert peptide sequences to masses, read from a Fasta file')
    parser.add_argument("filename", help="input Fasta file")
    parser.add_argument("--mass_type", choices=masstypes, default="a",
                        help="type of mass calculation: 'a' for average mass, 'm' for monoisotopic mass (default: 'a')")
    parser.add_argument("-c", "--charge", default="1", help="charge (default: 1)")
    args = parser.parse_args()

    # Extract arguments
    inputfile = args.filename
    masstype = args.mass_type
    charge = args.charge

    # Read sequences from input file
    order, seqs = fastaread(inputfile)

    # Set up output file name
    outputfile = f'pepmasses_from_{inputfile}'

    # Write peptide masses to output file
    with open(outputfile, 'w') as ofile:
        seq_number = 0
        for seq_name, seq in seqs.items():
            mass = pep2mass(seq, masstype)
            mass_to_charge = round((float(mass) + float(charge)) / float(charge), 4)
            # Write sequence name, peptide mass, charge, and sequence to output file
            ofile.write(f"{seq_name.split(',')[0]} {seq_name.split(',')[1]} "
                        f"{mass_to_charge} {charge} 0 {seq}\n")

    print("Peptide mass calculation completed. Results written to", outputfile)
