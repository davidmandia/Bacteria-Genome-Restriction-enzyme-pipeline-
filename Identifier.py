#!/usr/bin/env python3

"""
Script to identify and count peptides within a specified mass-to-charge (m/z) range in a given data file.

Usage:
    python peptide_analyzer.py <fileName> -s <start> -e <end> [-b <binsize>] [-w <window>] [-ss <stepsize>] [-ma <massaccuracy>]

Arguments:
    <fileName>: Name of the file containing peptide data. The format must contain the name of the sequence, m/z charge, and sequence.
    -s, --start: Lower bound of m/z range (required).
    -e, --end: Upper bound of m/z range (required).
    -b, --binsize: Binsize in m/z units for Mode 2 (optional).
    -w, --window: Window size in m/z units for Mode 3 (optional).
    -ss, --stepsize: Step size in m/z units for Mode 3 (optional).
    -ma, --massaccuracy: Mass accuracy of the instrument, specified in Dalton or ppm (optional).

Output:
    Mode 1: Prints the number of peptides within the specified m/z range.
    Mode 2: Generates a file containing uniquely identified proteins and their details. Plots a histogram of the mass-to-charge distribution.
    Mode 3: Generates a file containing uniquely identified proteins using the sliding window method. Plots a graph of the mass-to-charge distribution using sliding windows.
"""

import argparse

# Function to read data from the file and filter peptides by mass range
def readData(file="mass_to_charge file.txt", lower=1000, upper=1500):
    dataFile = open(file, "r")
    
    # Skip the first line as it's a header
    first_line = dataFile.readline()
    peps = {}  # Dictionary to store peptides, indexed by a unique identifier

    # Read and process each line in the file
    for line in dataFile:
        if not line.startswith('#') and not line.startswith("\n"):  # Skip comments and empty lines
            (protein, pepnum, mass) = line.split()[0:3]  # Extract protein, peptide number, and mass
            pepseq = line.split()[5]  # Extract peptide sequence

            lower = float(lower)
            upper = float(upper)
            mass = float(mass)

            # Skip peptides outside the specified mass range
            if (mass < lower or mass > upper):
                continue

            # Create a unique identifier for each peptide
            pepid = "_".join(list((protein, pepnum, pepseq)))
            peps[pepid] = mass

    dataFile.close()  # Close the file
    return peps

if __name__ == "__main__":
    # Function to check for input range errors
    def error_check(range0, range1):
        if range0 < 0 or range1 <= 0:
            raise Exception("Sorry, mass-to-charge ratio cannot be below zero")
        if range0 == range1:
            raise Exception("Sorry, the two bounds cannot be the same")
        if range0 > range1:
            raise Exception("Sorry, the lower bound cannot be higher than the upper bound")
        return range0, range1

    # Set up argument parser with descriptions and default values
    parser = argparse.ArgumentParser(description="The script can be run in different modes. For example, mode 1 returns the number of peptides in a given mass range. It calculates the total number of peptides given lower and upper bounds for m/z values. Mode 2 returns uniquely identified peptides given a range and a binsize through the plotting of a histogram. Mode 3, similarly to mode 2, returns uniquely identified proteins through the usage of sliding windows. Mode 3 requires the input of a range, a window size, and a stepsize.")

    parser.add_argument("fileName", help="name of file containing data. The format must contain the name of the sequence, m/z charge, and sequence")
    parser.add_argument("-s", "--start", help="lower bound of m/z range", default=1000.0, type=float, required=True)
    parser.add_argument("-e", "--end", help="upper bound of m/z range", default=1500.0, type=float, required=True)
    parser.add_argument("-b", "--binsize", type=float, default=False, help="binsize in m/z units for Mode 2", required=False)
    parser.add_argument("-w", "--window", type=float, default=False, help="windowsize in m/z units for Mode 3", required=False)
    parser.add_argument("-ss", "--stepsize", type=float, default=False, help="stepsize in m/z units for Mode 3", required=False)
    parser.add_argument("-ma", "--massaccuracy", default=0.1, help="please insert the mass accuracy of the instrument. The mass accuracy can be given in Dalton or ppm. If given in ppm, there must not be any whitespace between the number and ppm, for instance, 50ppm. Also, if given in ppm, the ppm is calculated from the average of the protein weights within the range given", required=False)
    args = parser.parse_args()

    file = args.fileName
    range0 = round(args.start, 4)
    range1 = round(args.end, 4)
    binsize = args.binsize
    range0, range1 = error_check(range0, range1)

    # Calculate the number of bins if binsize is given
    if binsize:
        numbins = int((range1 - range0) / binsize)

    # Assign window size and stepsize
    W = args.window
    s = args.stepsize

    # Determine the decimal precision of the stepsize
    len_s = len(str(s)) - 2

    # Determine which mode to run based on the input arguments
    def determine_mode(binsize, W, s):
        Mode = 1  # Default mode

        # Check for invalid binsize and stepsize values
        if binsize < 0 or s < 0:
            raise Exception('Please provide positive values for binsize and/or stepsize')

        if binsize > 0:
            if not s:
                Mode = 2
                print("You have picked Mode 2. It will return the identified proteins along with masses and ranges used to identify them")
            if s > 0:
                Mode = 3
                if not W:
                    print("Window size not specified. The program will use the binsize as window size")
                    W = binsize
                    Mode = 3

        if s > 0 and W > 0:
            Mode = 3

        if not W and s > 0:
            raise Exception("You must provide a window size to run the program in mode 3")

        print("The program will run in mode:", Mode)
        return Mode, W

    pepmass = readData(file, float(range0), float(range1))
    Mode, W = determine_mode(binsize, W, s)

    # Determine mass accuracy based on input
    def mass_acc(mass_acc):
        if "ppm" in str(mass_acc):
            try:
                ppm = int(mass_acc.split('ppm')[0])
                m_a = round(max([float(i) for i in pepmass.values()]) * ppm / 1000000, 4)
            except ValueError:
                raise Exception("Sorry, you must insert either a number (Dalton), or a number followed by ppm. Otherwise, you can omit the mass accuracy, the program will default it to 0.1 Dalton")
        else:
            try:
                m_a = round(float(mass_acc), 4)
            except ValueError:
                raise Exception("Sorry, you must insert either a number (Dalton), or a number followed by ppm")
        return m_a

    mass_accuracy = mass_acc(args.massaccuracy)

    # Determine the decimal precision of the mass accuracy
    len_m_a = len(str(mass_accuracy)) - 2

    # Check for errors specific to modes 2 and 3
    def error_checking_modes(Mode, binsize, mass_accuracy, s, W):
        if Mode == 2:
            if binsize < mass_accuracy:
                raise Exception('The value of binsize cannot be lower than the instrument accuracy of {} Dalton'.format(mass_accuracy))
        if Mode == 3:
            if W < mass_accuracy:
                raise Exception('The value of window size cannot be lower than the instrument accuracy. Defaulted at 0.1 Dalton')
            if s > W:
                raise Exception('The value of stepsize cannot be larger than the window')
            if s < 0.0001:
                raise Exception('The stepsize cannot be lower than the mass decimal point of 0.0001')

    error_checking_modes(Mode, binsize, mass_accuracy, s, W)

    # Execute the program based on the determined mode
    def execute_programme(Mode):
        if Mode == 1:
            print("There are", len(pepmass), "peptide within the selected range")

        import matplotlib.pyplot as plt
        if Mode == 2:
            output_file = 'identified_protein_from_inputfile_{}_in_mode_{}.txt'.format(file, Mode)
            outfile = open(output_file, 'w')
            numbins = int((range1 - range0) / binsize)
            bin_range_count = {"count": [], "label": []}
            prot_uniq = []

            for bin in range(0, numbins):
                lower = round(range0 + bin * binsize, len_m_a)
                upper = round(lower + binsize, len_m_a)
                count = 0

                for pep in pepmass.keys():
                    mz = float(pepmass[pep])
                    if (mz >= lower and mz < upper):
                        count += 1

                bin_range_count['count'].append(int(count))
                bin_range_count['label'].append((lower + upper) / 2)

                if count == 1:
                    for pep in pepmass.keys():
                        mz = float(pepmass[pep])
                        if (mz >= lower and mz < upper):
                            prot_uniq.append([pep, mz, lower, upper])

            print("Uniquely identified protein number:", len(prot_uniq))
            for uniq in prot_uniq:
                text_file_line = "The protein {prot} has been identified unique with mass {mass} within the range {lower} - {upper}\n".format(prot=uniq[0], mass=uniq[1], lower=uniq[2], upper=uniq[3])
                outfile.write(text_file_line)
            print('Uniquely identified proteins have been added to the file {}'.format(output_file))
            outfile.close()

            y = bin_range_count["count"]
            x = bin_range_count['label']

            plt.xlabel("mass to charge ratio (m/z)")
            plt.ylabel("frequency")
            plt.xlim([range0, range1])
            plt.plot(x, y)
            plt.savefig("Mass to charge distribution.png")
            plt.show()

        if Mode == 3:
            lower_slid = round(range0, len_s)
            upper_slid = round(lower_slid + W, len_s)
            numsteps = int((range1 - range0) / s)
            sliding_wind = {"bounds": [], "count": []}
            prot_uniq_sliding = {"prot": [], "mz": [], "lower": [], "upper": []}

            output_file = 'identified_protein_from_inputfile_{}_in_mode_{}.txt'.format(file, Mode)
            outfile = open(output_file, 'w')

            for step in range(0, numsteps):
                if upper_slid >= range1:
                    break

                count = 0
                for pep in pepmass.keys():
                    mz = float(pepmass[pep])
                    if (mz >= lower_slid and mz < upper_slid):
                        count += 1

                sliding_wind["count"].append(int(count))
                sliding_wind["bounds"].append((lower_slid + upper_slid) / 2)

                if count == 1:
                    for pep in pepmass.keys():
                        mz = float(pepmass[pep])
                        if (mz >= lower_slid and mz < upper_slid):
                            if pep not in prot_uniq_sliding["prot"]:
                                prot_uniq_sliding["prot"].append(pep)
                                prot_uniq_sliding["mz"].append(mz)
                                prot_uniq_sliding["lower"].append(lower_slid)
                                prot_uniq_sliding["upper"].append(upper_slid)

                lower_slid = round(lower_slid + s, len_s)
                upper_slid = round(lower_slid + W, len_s)

            print("Uniquely identified protein number:", len(prot_uniq_sliding["prot"]))
            for i in range(0, len(prot_uniq_sliding["prot"])):
                text_file_line = "The protein {prot} has been identified unique with mass {mass} within the range {lower} - {upper} in the sliding window method\n".format(prot=prot_uniq_sliding["prot"][i], mass=prot_uniq_sliding["mz"][i], lower=prot_uniq_sliding["lower"][i], upper=prot_uniq_sliding["upper"][i])
                outfile.write(text_file_line)
            print('Uniquely identified proteins have been added to the file {}'.format(output_file))
            outfile.close()

            y = sliding_wind["count"]
            x = sliding_wind["bounds"]

            plt.xlabel("mass to charge ratio (m/z)")
            plt.ylabel("frequency")
            plt.xlim([range0, range1])
            plt.plot(x, y)
            plt.savefig("Sliding window.png")
            plt.show()

    execute_programme(Mode)
