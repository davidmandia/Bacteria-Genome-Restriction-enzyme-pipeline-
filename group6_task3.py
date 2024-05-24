# This script written by Yaoran Yu, ID: 11349270, in group 6
# Define functions for reading files
import os

def fastaread(filename):
    with open( filename, "r") as fastafile:
        # Open the file and read all lines
        lines = fastafile.readlines();
    
    order = []
    seqs = {}
    for line in lines:
        if line.startswith('>'):
            seq_name = line.split()[0][1:] + ',' + line.split()[2]    # If the line starts with '>', get the sequence name
            order.append(seq_name)    # Split the name line to the name and the order number, use ','to distinguish
            seqs[seq_name] = ''    # Add the name and the order nunber to the list
        else:
            seqs[seq_name] = line.strip()    # Add the whole sequence line, which contains the sequence
    return (order, seqs)

# Define the function of converting peptide sequences into mass
def pep2mass(seq, type):
    mono = {'A' :71.0371, 'C' :103.0092, 'D' :115.0269, 'E' :129.0426, 'F' :147.0684, 'G' :57.0215, 'H' :137.0589, 'I' :113.0841, 'K' :128.0950, 'L' :113.0841, 'M' :131.0405, 'N' :114.0429, 'P' :97.0528, 'Q' :128.0586, 'R' :156.1011, 'S' :87.0320, 'T' :101.0477, 'V' :99.0684, 'W' :186.0793, 'Y' :163.0633, '*' :0.0}
    aver = {'A' :71.08, 'C' :103.14, 'D' :115.09, 'E' :129.12, 'F' :147.18, 'G' :57.05, 'H' :137.14, 'I' :113.16, 'K' :128.17, 'L' :113.16, 'M' :131.19, 'N' :114.10, 'P' :97.12, 'Q' :128.13, 'R' :156.19, 'S' :87.08, 'T' :101.10, 'V' :99.13, 'W' :186.21, 'Y' :163.18, '*' :0.0}
    waterM = 19.0106
    waterA = 19.0153
    
    if type == "mono":
        mass_dict = mono
        water = waterM
    else:
        mass_dict = aver
        water = waterA
    # Select mass dictionary and water mass based on type choice

    mass = round(sum(mass_dict[amino_acid] for amino_acid in seq) + water, 4)    # Calculate sequence quality to 4 decimal places

    return mass

if __name__ == "__main__":    # If run directly
    import argparse
    masstype = ['a','m']
    # Add parameters
    parser = argparse.ArgumentParser(description='convert peptide sequences to masses, read from a Fasta file')
    parser.add_argument("filename",help="The input filename, must be in fasta format", default="dummy.fasta")
    parser.add_argument("type",help="The masstype, choose a for average, m for monoisotopic", choices=masstype, default="a")
    parser.add_argument("output_name",help="The output filename", default="pepmasses.out")
    parser.add_argument("charge",help="The charge number", default="1")
    
    # Parsing parameters
    args = parser.parse_args()
    inputfile = args.filename
    masstype = args.type
    outputfile = args.output_name
    z = args.charge

    # Check if the input file exists
    if not os.path.isfile(inputfile):
        print(f"Error: The file '{inputfile}' does not exist.")
        exit(1)
    # Check if the masstype is valid
    if masstype not in ['a', 'm']:
        print(f"Error: The masstype '{masstype}' is not valid. Choose 'a' for average, 'm' for monoisotopic.")
        exit(1)
    # Check if the charge is a number
    try:
        float(z)
    except ValueError:
        print(f"Error: The charge '{z}' is not a number.")
        exit(1)
    order, seqs = fastaread(inputfile)
    # Open output file
    with open(outputfile,'w') as ofile:
        list = []
        list.append(['Prot_name','peptide','mass-to-charge','z','p','sequence'])
        # Subtitle
        for seq_name in seqs:
            seq = seqs[seq_name]    # Sequence
            mass = pep2mass(seq, masstype)    # Calculate the mass
            mass_to_charge = round((float(mass) + float(z))/float(z), 4)    # Calculate mass to charge
            list.append([seq_name.split(',')[0],seq_name.split(',')[1],mass_to_charge,z,'0',seq])    # Put results into a list, use ','to distinguish
        for row in list:    # Format the list and write into the file
            for item in row:
                text = str(item).center(20)    # Use center to separate the data, which is convenient for the next step to distinguish using spaces
                ofile.write(text)    # Write the data into the file
            ofile.write('\n')
else:
    print("run as module\n")    # If run as module
