#!/bin/bash

# Define variables
INPUT_FILE="genome.fasta"

min_orf=$1
restriction_type=$2
pepmin=$3
mass_type=$4
start=$5
end=$6

# Run script 1
python finding_orf.py -m "$min_orf" "$INPUT_FILE" > output1.txt

# Run script 2 (restr_enzyme.py)
python restr_enzyme.py -e "$restriction_type" -p "$pepmin" output1.txt > output2.txt

# Run script 3 (mass_Calculator.py)
python mass_Calculator.py --mass_type "$mass_type" output2.txt > output3.txt

# Run script 4 (Identifier.py)
python Identifier.py -s "$start" -e "$end" output3.txt > final_output.txt
echo "Final output file 'final_output.txt' has been generated."
echo "in mode 2 and 3 graph will be generated, but it has not been implemented "

# Clean up intermediate files
rm -f output1.txt output2.txt output3.txt
