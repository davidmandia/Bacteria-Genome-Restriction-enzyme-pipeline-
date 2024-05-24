



def readData(file="mass_to_charge file.txt", lower=1000, upper=1500):
    dataFile = open(file, "r")
    
    #skip the fist line as header
    fist_line = dataFile.readline()
    peps = {
    }  # dictionary containing the masses, indexed via a unique identifier

    ## assumes file format like:   YCL049C    5     774.4059  0 0 GVAMGNVK

    for line in dataFile:
        # I have added another condition to avoid error with empty lines
        if not line.startswith('#'):
            #remove the presence of new line character 
            if not line.startswith("\n"):
                # split the line, and assigned values of the list to different variables
                # Missed cleavage variables has been ignored 
                #print(line.split())

                (protein, pepnum, mass) = line.split()[0:3]
                pepseq = line.split()[5]

                #print("data = ",protein, pepnum, mass, pepseq)
                # print(mass)
                # print(lower)
                # print(type(mass), type(lower))

                # When run as pipelin in Bash it needs to be converted 
                lower = float(lower)
                upper = float(upper)
                mass = float(mass)
                print(type(mass))
                if (float(mass) < lower or float(mass) > upper):  # Exclude the proteins which mass-to charge ratio is not within the bounds
                    continue

#
# make a unique id for future reference for each peptide

                pepid = "_".join(list((protein, pepnum, pepseq)))
                peps[pepid] = mass
#Closed the file as per common practice 
    dataFile.close()
    return (peps)


if __name__ == "__main__":
    
    def error_check(range0, range1):
        
            #Raising errors with the input of data 
            ## Ranges cannot be negative

            # The lower bound can be 0, but not the upper bound 
            if range0 < 0 or range1 <= 0:
                raise Exception("Sorry,  mass-to charge ratio cannot be below zero")

            # Ranges cannot be the same 

            if range0 == range1:
                raise Exception("Sorry, They two bounds cannot be the same")

            # The upper bound cannot be lower than the lower bound

            if range0 > range1:
                raise Exception("Sorry, The lower bound cannot be higher than the upper bound")
            return range0, range1





    import argparse
    # The help message
    parser = argparse.ArgumentParser( description = "The script can be run in different modes. For example, mode 1 returns the number of peptides in a given mass range. It calculates  the total number of peptides given a lower and upper bounds for m/z values. Mode 2 return uniquely identified peptides given a range, and a binsize through the plotting of an histogram. Mode 3, similarly to mode 2, returns uniquely identified protein through the usage of sliding windows. Mode 3 requires the input of a range, a window size, and a stepsize.") 



    # Start and End arguments are mandatory. If no other arguments is given, the programm will default to Mode 1 

    parser.add_argument("fileName", help="name of file containing data. The format must contain in name of the sequence, m/z charge, and sequence")     # a positional argument
    parser.add_argument("-s", "--start",  help="lower bound of m/z range", default=1000.0, type=float, required=True)

    parser.add_argument("-e", "--end",    help="upper bound of m/z range", default=1500.0, type=float, required=True)



    #If binsize not given, it will be  defaulted to False. which will then switch to mode 1 



    parser.add_argument("-b", "--binsize", type= float, default = False ,help="binsize in m/z units for Mode 2", required=False) 

    #Both w and s are defaulted to False. if not given, to defualt to mode 1 or 2


    parser.add_argument("-w", "--window", type= float, default = False, help="windowsize in m/z units. This is the width of the window in mode 3", required=False)# Added the window size for mode 3 

    parser.add_argument("-ss", "--stepsize",  type= float , default = False , help="stepsize in m/z units. This is the size of the step in mode 3. The windonw will move to the right by stepsize Dalton", required=False)# Added the stepsize for mode 3 

    # Adding the mass accuracy. It is not required though. It can be given in Dalton or ppm 
    parser.add_argument("-ma", "--massaccuracy", default = 0.1 , help="please insert the mass accuracy of the instrument. The mass accuracy can be given in Dalton as number or in ppm. If given in ppm, There must not be any whitespace between the number an ppm, for instance, 50ppm. Also, if given in ppm, the ppm is calculated from the average of the protein weights within the range given ", required=False)
    args = parser.parse_args()


    file = args.fileName

    #the two ranges are approximated to the 4th decimal place. That is same decimal places of the mass of the peptides

    range0 = round(args.start, 4)
    range1 = round(args.end, 4)
    binsize = args.binsize
    range0, range1 = error_check(range0, range1)


    # the number of bins is calculated only if binsize is given


    if binsize == True: # Binsize is defaulted to False if not given. If given as number different from 0, binsize will be true

        #the number of bins is the range divided by the width of the bin
        # The number of bins must be an integer
        numbins = int((range1-range0)/binsize)

   
    # Window size and stepsize are assigned    
    W =  args.window
    s = args.stepsize

    # This will save the decimal place of the stepsize as number.
    
    
    #Useful for later approximations as we are working with floats additions 

    len_s = len(str(s)) - 2
    #print("Stepsize", len_s, type(len_s))



    def determine_mode(binsize, W, s):
        Mode = 1 # default mode

        # If given, stepsize or binsize cannot be negative 
        if binsize < 0 or s < 0:
            raise Exception('Please provide positive values for binsize and/or stepsize')

        # If the binsize is given and positive, and the stepsize not given. Mode 2 is used
        if binsize > 0 :
            if s == False:
                Mode = 2 
                print( "You have picked Mode 2. It will return  the identified proteins along with masses and ranges used to indetify them" )

        # The presence of stepsize will switch the programm to mode 3 , sliding window         
            if s > 0:
                Mode = 3

                # For usability the programme will use binsize as window if window size  is not given 
                if W == False:
                    print("Window size not specified. The programm will use the binsize as window size")
                    W = binsize
                    Mode  = 3


        # if stepsize and window size are given, the program will run in Mode 3 
        if s > 0 and W > 0:
            Mode =3






        if W == False and s > 0:
            raise Exception("you must provide a Window size to run the program in mode 3")






        print("The programm will run in mode:", Mode)
        
        # Added W to the return statement. Window can default to binsize, if window not given bus still on mode 3
        return Mode, W






    pepmass = {}

    #read file given file data and bounds
    pepmass = readData(file, float(range0) , float(range1))
    Mode, W = determine_mode(binsize, W, s)

    # mass _accuracy 
    # The user can provide the mass accuracy in ppm or Dalton
    def mass_acc(mass_acc):
        #The programme will recognise either a number or a number followed by ppm 

        if "ppm" in str(mass_acc):
            #This means the user gave the mass accuracy in parts per milion


            #Takes the fist part of the input and converts into int. This assumes the input is given as int followed by 'ppm'

            # print('mass in ppm =', ppm )
            try:

                ppm = int(mass_acc.split('ppm')[0])

                # The programme will use the highest weight within the bounds to calculate the mass accuracy

                # The highest value has been chosen as that will give the highest mass accuracy in Dalton. Which means an higher confidence in measurements


                m_a = round(max([float(i) for i in pepmass.values()]) * ppm / 1000000, 4)
            except ValueError:
                raise Exception("Sorry, You must insert either a number (Dalton), or a number followed by ppm. Otherwise, you can omit the mass accuracy, the program will default it to 0.1 Dalton")  

        else:
            try:
                m_a = round(float(mass_acc),4)
            except ValueError:
                raise Exception("Sorry, You must insert either a number (Dalton), or a number followed by ppm") 
            #print('Mass in Dalton', m_a)

        return m_a
       

    # mass accuracy for mode 2 and 3 
    mass_accuracy = mass_acc(args.massaccuracy)

    # This will used to round up to the same decimal place as the mass accuracy
    
    len_m_a = len(str(mass_accuracy)) - 2


    def error_checking_modes(Mode, binsize, mass_accuracy, s, W):

        if Mode == 2:
            # Bin size cannot be lower than the instrument precision 
            if binsize < mass_accuracy :
                raise Exception('The values of binsize cannot be lower than the instrument accuracy of {} Dalton'.format(mass_accuracy))

        if Mode == 3:
            #Window size cannot be lower than the instrument precision 
            if W < mass_accuracy :

                raise Exception('The values of window size cannot be lower than the instrument accuracy. Defaulted at 0.1 Dalton')


            if s > W:
                raise Exception('The value of stepsize cannot be larger than the window')

                # Stepsize must be higher than the lowest decimal place of the m/z ratio
            if s < 0.0001:
                raise Exception('The stepsize cannot be lower than the mass decimal point of 0.0001')
            
            
            
    # Function called without being assigned to a variable as it is just checking values for that specific mode
            
    error_checking_modes(Mode, binsize, mass_accuracy, s, W)      
    

    #Plain Mode 1 
    
    
    def execute_programme(Mode):
        
        if Mode == 1:
             #Return the number of peptides with that m/z range 
            print("There are", len(pepmass), "peptide within the selected range")
        # For mode 2 and 3 an histogram is drawn
        
        import matplotlib.pyplot as plt
        if Mode == 2:
            #Output file for uniquely identified protein
            output_file ='identified_protein_from_inputfile_{}_in_mode_{}.txt'.format(file, Mode)
            outfile= open(output_file, 'w')
            numbins = int((range1 - range0)/binsize)
            bin_range_count = {
            "count": [],
            "label": []
                }
            # Create an empty multi dimensional matrix which will contain the peptide name, mass, and the bounds which made uniquely identified 
            prot_uniq = []

    # Itirate through the bins 




            for bin in range(0,numbins):

            # Upper and lower bounds for that specific bin
                lower = round(range0 + bin*binsize, len_m_a)

                # The round function is used to avoid error when summing between floats
                upper = round(lower + binsize, len_m_a)


                count = 0 
                # Itirate through the pepmass dictionary 
                for pep in pepmass.keys():
                    #Obtain the mass to charge value for that specific pep
                    mz = float(pepmass[pep])

                    # See if the m/z is within the boundaries of that specific bin 

                    if ( mz >= lower and mz < upper):
                        count += 1

                # Append the values of the count to the dictionary for later use as y  

                bin_range_count['count'].append(int(count))

                # Append the average of the bin bounds to the dictionary to be later used on the x axis 
                bin_range_count['label'].append((lower+upper )/2)
                # The programme will return the identified protein along with masses and ranges  

            # Find the peptide that are unique to that one bin and append the to the multidimensional matrix 
                if count == 1:
                    for pep in pepmass.keys():
                        mz = float(pepmass[pep])
                        if ( mz >= lower and mz < upper):
                            prot_uniq.append([pep, mz, lower, upper])


            # Itirate through the uniq prot matrix
            print("Uniquely identified protein number:", len(prot_uniq))                
            for uniq in prot_uniq:
                text_file_line = "The protein {prot} has been identified unique with mass {mass} within the range {lower} - {upper}\n".format(prot=uniq[0], mass=uniq[1], lower=uniq[2], upper=uniq[3])
                
                
                
                outfile.write(text_file_line)
            print('Uniquely identified proteins have ben added to the file {}'.format(output_file))
            outfile.close()
                
            # This part is for the histogram

            y = bin_range_count["count"]

            x = bin_range_count['label']


            # Labels

            plt.xlabel("mass to charge ratio ( m/z)")
            plt.ylabel("Frequency")

            # limits on x
            plt.xlim([range0, range1])

            plt.plot(x, y)
            plt.show()
            plt.savefig("Mass to charge distribution.png")

        # Mode 3  - sliding window 
    

        if Mode == 3:
            
            
            output_file ='identified_protein_from_inputfile_{}_in_mode_{}.txt'.format(file, Mode)
            outfile= open(output_file, 'w')

            # Sliding window data for the graph 

            sliding_wind = {
                "bounds": [],
                "count": []
            }

            num_steps = int((range1 -range0) / s)

            # The first window is exacty the bin 



            lower_slid = round(range0 , len_s )
            upper_slid = round(lower_slid + W, len_s )


            prot_uniq_sliding = {
                "prot": [],
                "mz": [],
                "lower": [],
                "upper": []
            } # Empty matrix for protein identified via sliding window 


            # Itirate through the number of steps within the range 
            for i in range(0,num_steps):



                # Upper and lower bounds for that specific bin    



                if upper_slid >= range1:
                    break

                count = 0 # Count is re-started to 0 for every window
                # Itirate through the pepmass dictionary 
                for pep in pepmass.keys():
                    #Obtain the mass to charge value for that specific pep
                    mz = float(pepmass[pep])

                    # See if the m/z is within the boundaries of that specific bin 

                    if ( mz >= lower_slid and mz < upper_slid):
                        count += 1

                # Append the values of the count to the dictionary for later use as y  

                sliding_wind['count'].append(count)

                # Append the average of the bin bounds to the dictionary to be later used on the x axis 
                sliding_wind['bounds'].append((lower_slid+upper_slid )/2)


                if count == 1:
                    for pep in pepmass.keys():
                        mz = float(pepmass[pep])
                        if ( mz >= lower_slid and mz < upper_slid and pep not in prot_uniq_sliding["prot"]):
                            prot_uniq_sliding["prot"].append(pep)
                            prot_uniq_sliding["mz"].append(mz)
                            prot_uniq_sliding["lower"].append(lower_slid)
                            prot_uniq_sliding["upper"].append(upper_slid)

                lower_slid = round(lower_slid + s, len_s)
                upper_slid = round(lower_slid + W, len_s)
                #print("lower-slid", lower_slid,"upper_slid", upper_slid)



            # Itirate through the uniq prot matrix
            print("Uniquely identified protein number:", len(prot_uniq_sliding['prot']))                
            for i in range(len(prot_uniq_sliding["prot"])):
                #Write output for Mode 3 
                text_file_line = "The protein {prot} has been identified unique with mass {mass} within the range {lower} - {upper} in the sliding window method\n".format(prot=prot_uniq_sliding["prot"][i], mass=prot_uniq_sliding["mz"][1], lower=prot_uniq_sliding["lower"][i], upper=prot_uniq_sliding["upper"][i])
                
                
                
                outfile.write(text_file_line)
            print('Uniquely identified proteins have ben added to the file {}'.format(output_file))
            outfile.close()



            y_slid = sliding_wind["count"]

            x_slid = sliding_wind['bounds']


            # Labels

            plt.xlabel("mass to charge ratio ( m/z) sliding window ")
            plt.ylabel("Frequency")

            # limits on x
            plt.xlim([range0, range1])

            plt.plot(x_slid, y_slid)

            plt.savefig("sliding_window.png")
            plt.show()
    
    #Execute the programme
    execute_programme(Mode)

else:
    print("run as module\n")



        
        
