#!/usr/bin/env python
#
# protein digestion tool
#
#####  Mike Frost 11387111  Task 2 Group 6
#

import re, sys, os  # not essential for this version, but you might need them

#---------------------------------------------------
#
#        Fasta file read function
#
def fastaread(filename):
#
    """Opens and reads in a Fasta formatted file 
    """
    try:
        with open( filename, "r") as fastafile: # The 'with' operation allows us to auto-close the file
            lines = fastafile.readlines(); # lines is declared as a list containing the lines in the file
    except:
        print("File could not be opened. Check permissions")


#

    order = []    	# declare a list for sequence names
    seqs = {}		# declare a dictionary of sequences
    seqs_temp = ""			 # temporary string to build sequence
    seqs_name_temp1 = ""	 # temporary string to build sequence name
    seqs_name_temp2 = ""	 # temporary string to build sequence name
    seqs_name = ""			 # string to hold final sequence name for the order list

    for line in lines:     			# iterate over fasta file, line by line
        
        if line.startswith('>'):	# for header lines which start with '>'
            seqs_temp = "" 			# Set temporary sequence string variable to empty because we are starting a new sequence
            seq_name_temp1 = line.split(" ")[0]		# Drop any characters to the right of the first space from the sequence name.
            if ',' in seq_name_temp1:				# Check for a comma in the first part of the sequence name. Seen in some early test data
                seq_name_temp2 = seq_name_temp1.split(",")[1]	# If found, drop any characters to the left of the comma
                seq_name = '>' + seq_name_temp2             	# Reinstate the fasta header symbol removed by the previous line at the beginning of the string
                order.append(seq_name)   						# Add the sequence name to the order list variable       
            else: # No comma so no adjustment to string
                seq_name = seq_name_temp1		# Adopt the sequence name string without further amendment
                order.append(seq_name)			# Add the sequence name to the order list variable

        else:		# If not a header line must be an amino acid sequence line.
            if line.rstrip().endswith('*'):			# Lines ending with * (stop) are the last line. If so do the following:
                seqs_temp = seqs_temp + line.rstrip()	# Strip any new line characters from sequence
                seqs[seq_name] = seqs_temp		# Add the sequence to the dictionary, keyed on it's name
            else: # It must be a continuing sequence that wraps
                seqs_temp = seqs_temp + line.rstrip()	# Append to the sequence. Still strip new line characters 

#    print ('%d sequences read in' % len(order))  # Uncomment to DEBUG 

    return (order, seqs)

# -------------------------------------------------
#
#     Digest function
#

def digest(sequence, enzyme):

    peps = []   # A list to contain the peptides generated from a sequence
    pep_temp = "" 			 # temporary string to build peptide sequence
    seq_ct = 0 			 # initialise a loop iterator for processing the sequence
    if enzyme == 't':  			 # If the enzyme is trypsin
        cut_list = ['K','R']  	 # cut after lysine or arginine
    elif enzyme == 'r':  		 # If the enzyme is endoproteinase Arg-C     
        cut_list = ['R']  		 # cut after arginine
    elif enzyme == 'k':  		 # If the enzyme is endoproteinase Lys-C       
        cut_list = ['K']  		 # cut after lysine 
    elif enzyme == 'e':  		 # If the enzyme is V8 proteinase (Glu-C)       
        cut_list = ['E']  		 # cut after glutamic acid
                    
    while seq_ct < len(sequence):			# loop until we have reached the end of the sequence
        if sequence[seq_ct] == '*':			# If we have reach the 'stop' then the following actions
            peps.append(pep_temp)			# append the last peptide to the peps list variable
            pep_temp = ""					# empty the temporary string for safety
        pep_temp = pep_temp + sequence[seq_ct]	# Add the current amino acid to the temporary peptide sequence
        if sequence[seq_ct] in cut_list: 		# Check to see if it's a cut point
            if sequence[seq_ct+1] != 'P':		# Check to make sure the cut is not prevented by a proline residue in the next position ie. if not proline
                peps.append(pep_temp) 			# Not proline so the cut occurs, append the temporary sequence to peptide list
                pep_temp = ""					# empty the temporary string (essential here)
        seq_ct += 1 							# End of loop so increment count by one
        
#    print (peps)  # Uncomment to DEBUG
    return (peps)  		# Return list of peptides

# -------------------------------------------------
#
#     Missed cleavages function
#

def mclpeps(peptide, pep_ct, miss_cl): # define the missed cleavages function as the following
    pep_mc = []       # Initialise a list variable to return from the function
    pep_mc_temp = ""  # Initialise a temporary peptide variable
    miss_cl_ct = 0   	# declare a missed cleavage loop count variable
    while miss_cl_ct <= miss_cl and miss_cl_ct + pep_ct < len(peps): # loop while missed cleavage count is less than the limit and we have not gone past the end of the list
        if miss_cl_ct == 0:        # check if it is the first iteration
            pep_mc_temp = peptide  # on first iteration populate temporary peptide string with peptide
        else:                      # on subsequent iterations
            pep_mc_temp = pep_mc_temp + peps[miss_cl_ct + pep_ct]  # add next peptide
        miss_cl_ct += 1  # At the end of the while loop add one to the missed cleavage count
        pep_mc.append(pep_mc_temp) # append the 'missed cleavages' peptide to the return variable
    miss_cl_ct = 0      # Clear the missed cleavage count ready for the next peptide 
#    print(pep_mc)
    return (pep_mc) # return list created by the function to main program


#-----------------------------------------------------
#
# main program
#

if __name__ == "__main__":
 
#
#  We've given you some argparse code to get you started, but you can adapt this and add more if you wish
#

    import argparse # import the argument parser library

    enzymes = {"t" : "trypsin",
               "r" : "endoproteinase Arg-C",
               "k" : "endoproteinase Lys_C",
               "e" : "V8 proteinase (Glu-C)"
              } # enzyme dictionary
    
    pepmins = []   #  create pepmins as a list and populate it with a loop
    it = 1
    max_it = 10  
    while it <= max_it:        #  valid peptide lengths from 1 to 10.
        pepmins.append(it)
        it += 1

    missed = [0, 1, 2]   #  create a list of valid choices for missed cleavages

    parser = argparse.ArgumentParser(description='digest proteins read in from a Fasta file', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("infileName")     # first positional argument input file name
    parser.add_argument("-of", "--outfileName", type=str, default='to_screen', help='''Please specify an output file name. 
Output goes to screen if none is supplied''', metavar='<str>')  # optional argument output file name. Will output to screen if not supplied
    parser.add_argument("-e", "--enzyme", help='''specify the enzyme from one of
t trypsin,
r endoproteinase Arg-C,
k endoproteinase Lys-C,
e V8 proteinase (Glu - C))''', choices=enzymes.keys(), default='t')   # optional argument to specify enzyme. Default is trypsin
    parser.add_argument("-p", "--pepmin", default=5, type=int, help='''specify a minimum peptide length from 1 to 10 to
include in the results (default: 5)''', choices=pepmins, metavar='<int>')  # optional argument to specify minimum peptide length. Default is 5
    parser.add_argument("-m", "--missed", default=0, type=int, help='''specify a maximum number of missed cleavages from
0 to 2 to allow for simple simulation of
partial digestion (default: 0)''', choices=missed, metavar='<int>')  # optional argument to specify maximum missed cleavages. Default is 0
    
    args = parser.parse_args()

    inputfile = args.infileName    # assign input file argument to variable
    enzyme = args.enzyme           # assign enzyme argument to variable
    print("The digestion enzyme is", enzymes[enzyme]) # Message for screen only
    pepmin = args.pepmin           # assign minimum peptide length argument to variable
    print("The minimum peptide length in the results is", pepmin, "amino acids")  # Message for screen only
    miss_cl = args.missed          # assign missed cleavages argument to a variable
    print("All permutations up to", miss_cl, "missed cleavage(s) will be included")  # Message for screen only
    outputfile = args.outfileName  # assign output file argument to variable
    if outputfile =='to_screen':
        print("No output filename given. Output will be to screen only\n")
    else:
        print("Output will be to file", outputfile)
        ofile = open(outputfile,'w')   # open a file handle for writing the output file
    
#
    order = []    	# declare order list in main
    seqs = {}		# declare seqs dictionary in main
    if os.path.isfile(inputfile): # if the input file exists
        (order, seqs) = fastaread(inputfile)   # process the input file
    else: # raise the missing file issue
        print("The input file", inputfile, "does not exist. Check file name.")
    peps = []   	# declare peptides list in main
    peps_out_ct = 1  	# declare a peptide output count variable to use in the output
    pep_mc_temp = "" # declare a temporary string variable to hold the peptide and grow when missed cleavages are set to more than 0. 

    for prot in order:  # loop through the list of proteins from the input file
#        print(order) # Uncomment for DEBUG data to screen
        peps = digest(seqs[prot], enzyme)  # call the digest function for each protein
#        print(peps) # Uncomment for DEBUG data to screen
        for pep_ct, peptide in enumerate(peps):     # loop through the peptides produced for each protein produced by the digest function
            pep_mc = mclpeps(peptide,pep_ct,miss_cl) # for each peptide call the missed cleavages peptides function to get all peptide combinations
            for frag_ct, pep_mc_frag in enumerate(pep_mc): # loop through the peptides from the missed cleavages function
                if len(pep_mc_frag) >= pepmin:  # if the peptide is longer than or equal to the minimum length do the following steps
                    if outputfile != 'to_screen': # If output is not to screen (default) write to file
                        ofile.write(f"{prot} peptide {peps_out_ct:03d} missed={frag_ct}\n{pep_mc_frag}\n")  # output the protein name and peptide to the output file
                        peps_out_ct += 1    # increment the petides output count
                    else:  # output to screen only
                        print(f"{prot} peptide {peps_out_ct:03d} missed={frag_ct}\n{pep_mc_frag}")  # output the protein name and peptide to standard output
                        peps_out_ct += 1    # increment the petides output count
        peps_out_ct = 1     # Reset the peptide output count ready for the next peptide
            
    if outputfile != 'to_screen': # if the output is to file
        ofile.close()      # close the output file handle

else:
    print("run as module\n")
