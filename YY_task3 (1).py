def fastaread(filename):
    with open( filename, "r") as fastafile:
        lines = fastafile.readlines();
    
    order = []
    seqs = {}
    for line in lines:
        if line.startswith('>'):
            seq_name = line.split()[0][1:] + ',' + line.split()[2]
            order.append(seq_name)
            seqs[seq_name] = ''
        else:
            seqs[seq_name] = line.strip()
    return (order, seqs)

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

    mass = round(sum(mass_dict[amino_acid] for amino_acid in seq) + water, 4)

    return mass

if __name__ == "__main__":
    import argparse
    masstype = ['a','m']
    parser = argparse.ArgumentParser(description='convert peptide sequences to masses, read from a Fasta file')
    parser.add_argument("filename", default="dummy.fasta")
    parser.add_argument("type", choices=masstype, default="a")
    parser.add_argument("output_name", default="pepmasses.out")
    parser.add_argument("charge", default="1")
        
    args = parser.parse_args()
    inputfile = args.filename
    masstype = args.type
    outputfile = args.output_name
    z = args.charge
    order, seqs = fastaread(inputfile)
    with open(outputfile,'w') as ofile:
        seq_number = 0
        list = []
        list.append(['Prot_name','peptide','mass-to-charge','z','p','sequence'])
        for seq_name in seqs:
            seq = seqs[seq_name]
            mass = pep2mass(seq, masstype)
            mass_to_charge = round((float(mass) + float(z))/float(z), 4)
            list.append([seq_name.split(',')[0],seq_name.split(',')[1],mass_to_charge,z,'0',seq])
        for row in list:
            for item in row:
                text = str(item).center(20)
                ofile.write(text)
            ofile.write('\n')
else:
    print("run as module\n")