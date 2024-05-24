#this takes a file and reads each line. If the line starts with >, it is stored as a header. If not, strips the line of whitespace (\n) and stores it in the line_collect variable. It then converts this to a string. 


def open_fastafile(file_name):
    line_collect = "" #empty strings
    header = ""
    with open(file_name, 'r') as fileObj:
        line = fileObj.readline() #reads lines one at a time
        for line in fileObj:
            if line.startswith(">"):
                header = line  # Store the line as header if begin with >
            else:
                line = line.strip()
                line_collect = line_collect + line.strip() #store stripped lines in line_collect
                genome_string = ''.join(line_collect) #store as string in new variable
    fileObj.close()
    return genome_string

gencode_dictionary = { "aaa": "K", "aac": "N", "aag": "K", "aat": "N", "aca": "T", "acc": "T",
"acg": "T", "act": "T", "aga": "R", "agc": "S", "agg": "R", "agt": "S","ata": "I", "atc": "I", "atg": "M", "att": "I", "caa": "Q", "cac": "H","cag": "Q", "cat": "H", "cca": "P", "ccc": "P", "ccg": "P", "cct": "P","cga": "R", "cgc": "R", "cgg": "R", "cgt": "R", "cta": "L", "ctc": "L","ctg": "L", "ctt": "L", "gaa": "E", "gac": "D", "gag": "E", "gat": "D","gca": "A", "gcc": "A", "gcg": "A", "gct": "A", "gga": "G", "ggc": "G","ggg": "G", "ggt": "G", "gta": "V", "gtc": "V", "gtg": "V", "gtt": "V","taa": "*", "tac": "Y", "tag": "*", "tat": "Y", "tca": "S", "tcc": "S","tcg": "S", "tct": "S", "tga": "*", "tgc": "C", "tgg": "W", "tgt": "C", "tta": "L", "ttc": "F", "ttg": "L", "ttt": "F"}

complement_dictionary = {
    "c": "g",
    "g": "c",
    "a": "t",
    "t": "a"
}

def revcomplement(sequence):
    reverse_sequence = sequence[::-1]  # Reverse the sequence
    sequence_str = "".join(str(reverse_sequence))  # Convert sequence to string
    sequences = ""
    reverse_complement = sequence_str.translate(str.maketrans(complement_dictionary))
    sequences += reverse_complement + "\n"  # Add the reverse complemented sequence to sequences
    return sequences

#overall: iterates over DNA in increments of 3 looking for a start codon starting at position frameoffset. If a start is found, sets start position and records the codon and matches them to their corresponding residues in the dictionary. If the orf is not smaller than the specified min length, it appends its residues, name, length and start position to lists.


import re

def find_orf(sequence, frameoffset, min_orf_length):
    start_codon = "atg"
    stop_codons = ["taa", "tag", "tga"]
    orf_residues = [] #list
    current_orf = "" #string
    start_index = []
    in_orf = False
    protein_counter = 1
    names = []
    lenghts = []
    starts = []
    orfs = []

    for i in range(frameoffset, len(sequence) - 2, 3): #iterates over dna starting at frameoffset position incrementing i by 3. It stops at -2 to stop 2 positions from the end sequence, ensuring the last codon is counted. 
        codon = sequence[i:i + 3] #specifies a codon is from position i to i+3
        start_index.append(i) #append start position of codon as i
        residue = gencode_dictionary.get(codon, "?") #looks for appropriate codon in the gencode_dictionary, if none found it adds a "?" using get()                                 
        if codon == start_codon and not in_orf: #sees if codon is a start and if already in an ORF. if both conditions are true, an orf starts
            in_orf = True 
            start_position= i #sets start position
            current_orf = residue #stores residues corresponding to codon 
            protein_name = f"CLAUD_transcript_00{protein_counter}" #f-string of name including counter 
        elif codon in stop_codons and in_orf: #checks for stop codon, if both true marks end of orf
            in_orf = False #takes you out of orf
            if len(current_orf) >= min_orf_length: #ensures captured orfs are not below minimum specified length
                end_position = i + 3 #counts end position 
                orf_length = end_position - start_position #calculates orf length 
                
                # The whole ORf is ignored if ? is present
                if "?" in current_orf:
                    print("There is an N nucleotide in the sequence. The programme will ignore the Open Reading Frame")
                else: #appends everything else to corresponding list if > min length and no n in ORF                  
                   
                    orf_residues.append(protein_name + str(start_position) + str(orf_length) + current_orf) #appends string containing start pos, length and orf residues to orf_residues 
                    orfs.append(current_orf)
                    names.append(protein_name)
                    lenghts.append(orf_length)
                    starts.append(start_position)
                    protein_counter += 1
            current_orf = "" #tells program to reset the current_orf string ready to store next orf
        elif in_orf:
            current_orf += residue #residue appeneded to current_orf, += concatenates the residues together as a string

    return [orf_residues, names, lenghts, starts, orfs]



if __name__ == "__main__":
    import argparse
    Frameoffset= [1, 2, 3]
    parser = argparse.ArgumentParser(description='find ORFs in a nucleotide genome sequence')
    parser.add_argument("fasta_file")      #three arguments specified, .fa file, min_orf_length and frameoffset
    parser.add_argument("-m", "--min_orf_length", help="the minimum length of an ORF in amino acids", type=int)
    parser.add_argument("-f", "--frameoffset", help="the starting position (1, 2 or 3) the ORFs are read from", type = int,choices = Frameoffset)
    args = parser.parse_args()
    inputfile = args.fasta_file 
    output_file = (f'ORF.fasta_frame{args.frameoffset}') #saves output to file
    outfile= open(output_file, 'w')
    genome=[]
    genome= open_fastafile(inputfile)
    if args.min_orf_length == None: #errors displayed and program terminated if args not added correctly 
        print("Error: specify the minimum length (integer) of open reading frame you want to search for after -m")
    if args.frameoffset == None:
        print("Error: specify the reading frame, either 1, 2 or 3, you to search for open reading frames after -f")
    orf_residues, names, lenghts, starts, orfs = find_orf(genome , args.frameoffset , args.min_orf_length)
    for i in range(len(orf_residues)): #writes details of orfs to file with specified format
         outfile.write(">{}_forward \t (Start: {} \t Length: {} bases) \n{}*\n".format(names[i],  starts[i], starts[i], orfs[i])) #{} are placeholders, first one holds names, second holds start pos and so on using the format(). It is a dynamic string
    orf_residues, names, lenghts, starts, orfs = find_orf(revcomplement(genome) , args.frameoffset , args.min_orf_length)
    for i in range(len(orf_residues)): #writes details of orfs to file with specified format
         outfile.write(">{}_reverse \t (Start: {} \t Length: {} bases) \n{}*\n".format(names[i],  starts[i], starts[i], orfs[i]))
            
            
outfile.close()


