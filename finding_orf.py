#!/usr/bin/env python3
"""
Script to find Open Reading Frames (ORFs) in a bacterial genome from FASTA file.

Usage:
    python find_orfs.py <fasta_file> -m <min_orf_length> 

Arguments:
    <fasta_file>: Input FASTA file containing the bacterial genome.
    -m, --min_orf_length: Minimum length of an ORF in amino acids.

Output:
    Generates a  file containing the found ORFs and their details.
"""

import argparse

# Genetic code dictionary
GENETIC_CODE = {
    "aaa": "K", "aac": "N", "aag": "K", "aat": "N",
    "aca": "T", "acc": "T", "acg": "T", "act": "T",
    "aga": "R", "agc": "S", "agg": "R", "agt": "S",
    "ata": "I", "atc": "I", "atg": "M", "att": "I",
    "caa": "Q", "cac": "H", "cag": "Q", "cat": "H",
    "cca": "P", "ccc": "P", "ccg": "P", "cct": "P",
    "cga": "R", "cgc": "R", "cgg": "R", "cgt": "R",
    "cta": "L", "ctc": "L", "ctg": "L", "ctt": "L",
    "gaa": "E", "gac": "D", "gag": "E", "gat": "D",
    "gca": "A", "gcc": "A", "gcg": "A", "gct": "A",
    "gga": "G", "ggc": "G", "ggg": "G", "ggt": "G",
    "gta": "V", "gtc": "V", "gtg": "V", "gtt": "V",
    "taa": "*", "tac": "Y", "tag": "*", "tat": "Y",
    "tca": "S", "tcc": "S", "tcg": "S", "tct": "S",
    "tga": "*", "tgc": "C", "tgg": "W", "tgt": "C",
    "tta": "L", "ttc": "F", "ttg": "L", "ttt": "F"
}

# Complement bases dictionary
COMPLEMENT_DICT = {
    "c": "g", "g": "c", "a": "t", "t": "a"
}

def open_fastafile(file_name):
    """
    Opens and reads the FASTA file.

    Parameters:
        file_name (str): The name of the FASTA file.

    Returns:
        str: The genome sequence read from the file.
    """
    genome_string = ""
    with open(file_name, 'r') as file_obj:
        for line in file_obj:
            if not line.startswith(">"):
                genome_string += line.strip()
    return genome_string

def revcomplement(sequence):
    """
    Generates the reverse complement of a DNA sequence.

    Parameters:
        sequence (str): The DNA sequence.

    Returns:
        str: The reverse complement of the input sequence.
    """
    reverse_sequence = sequence[::-1]
    reverse_complement = reverse_sequence.translate(str.maketrans(COMPLEMENT_DICT))
    return reverse_complement

def find_orf(sequence, min_orf_length):
    """
    Finds Open Reading Frames (ORFs) in a DNA sequence.

    Parameters:
        sequence (str): The DNA sequence.
        min_orf_length (int): The minimum length of an ORF in amino acids.

    Returns:
        tuple: Lists of ORF details including residues, names, lengths, and start positions.
    """
    start_codon = "atg"
    stop_codons = ["taa", "tag", "tga"]
    orf_details = []
    current_orf = ""
    in_orf = False
    protein_counter = 1
    names = []
    lengths = []
    starts = []

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        residue = GENETIC_CODE.get(codon, "?")
        if codon == start_codon and not in_orf:
            in_orf = True
            start_position = i
            current_orf = residue
            protein_name = f"CLAUD_transcript_00{protein_counter}"
        elif codon in stop_codons and in_orf:
            in_orf = False
            if len(current_orf) >= min_orf_length:
                end_position = i + 3
                orf_length = end_position - start_position
                if "?" not in current_orf:
                    orf_details.append((protein_name, start_position, orf_length, current_orf))
                    names.append(protein_name)
                    lengths.append(orf_length)
                    starts.append(start_position)
                    protein_counter += 1
            current_orf = ""
        elif in_orf:
            current_orf += residue

    return orf_details, names, lengths, starts

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find ORFs in a nucleotide genome sequence')
    parser.add_argument("fasta_file", help="input FASTA file")
    parser.add_argument("-m", "--min_orf_length", help="the minimum length of an ORF in amino acids", type=int)
    args = parser.parse_args()

    if args.min_orf_length is None:
        parser.error("Please specify the minimum length of an ORF using -m")

    output_file = 'ORF.fasta.txt'
    with open(output_file, 'w') as outfile:
        genome = open_fastafile(args.fasta_file)
        orf_details, _, lengths, starts = find_orf(genome, args.min_orf_length)
        for i, (protein_name, start_position, orf_length, current_orf) in enumerate(orf_details):
            outfile.write(">{}_forward \t (Start: {} \t Length: {} bases) \n{}*\n".format(protein_name, starts[i], lengths[i], current_orf))

        rev_comp_genome = revcomplement(genome)
        orf_details_rev, _, lengths_rev, starts_rev = find_orf(rev_comp_genome, args.min_orf_length)
        for i, (protein_name, start_position, orf_length, current_orf) in enumerate(orf_details_rev):
            outfile.write(">{}_reverse \t (Start: {} \t Length: {} bases) \n{}*\n".format(protein_name, starts_rev[i], lengths_rev[i], current_orf))
            print(">{}_reverse \t (Start: {} \t Length: {} bases) \n{}*".format(protein_name, starts_rev[i], lengths_rev[i], current_orf))