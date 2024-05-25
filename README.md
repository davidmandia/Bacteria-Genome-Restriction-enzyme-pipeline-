# Peptide Analysis Pipeline

This pipeline automates the process of identifying and analyzing peptides from a genome FASTA file. It involves four main steps:

1. **Finding ORFs**: Identifying open reading frames (ORFs) from the genome sequence.
2. **Digesting Sequences**: Using restriction enzymes to digest sequences into peptides.
3. **Calculating Peptide Masses**: Computing the masses of the digested peptides.
4. **Identifying Peptides**: Counting peptides within a specified mass-to-charge (m/z) range.

## Usage

### Step-by-Step Guide

1. **Finding ORFs**:
   - Run the `finding_orf.py` script to identify ORFs.
     ```bash
     python finding_orf.py -m 100 genome.fasta > output1.txt
     ```

2. **Digesting Sequences**:
   - Use the `restr_enzyme.py` script to digest sequences.
     ```bash
     python restr_enzyme.py -e t -p 5 output1.txt > output2.txt
     ```

3. **Calculating Peptide Masses**:
   - Run the `mass_Calculator.py` script to calculate peptide masses.
     ```bash
     python mass_Calculator.py --mass_type a output2.txt > output3.txt
     ```

4. **Identifying Peptides**:
   - Finally, use the `Identifier.py` script to identify peptides within the specified m/z range.
     ```bash
     python Identifier.py -s 1000 -e 1500 output3.txt > final_output.txt
     ```

### Unified Bash Script

To automate the process, use the provided bash script `run_pipeline.sh`:
```bash
./run_pipeline.sh <min_orf> <restriction_type> <pepmin> <mass_type> <start> <end>

### Example of output in Mode 2

   - It returns a graph with count of peptides fragment with m/z ratio within that bin
![Mass to charge distribution](https://github.com/davidmandia/Bacteria-Genome-Restriction-enzyme-pipeline-/blob/main/Mass%20to%20charge%20distribution.png)

