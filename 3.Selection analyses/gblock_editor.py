#! /usr/bin/python3

# edit gblock output files to fasta format

import MultiFastaLib as mfl
import glob

gbfiles = glob.glob("*gb")  #OG0006926_cds.aln-gb 
for i in gbfiles:
    newfile = i.replace("-gb","-gb.fas")
    seq_list = []
    with open(i, "r") as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                if ">" in line:
                    # new = line[0:5]

                    new = line # modified for non-single-copy orthologs

                else:
                    new = line.replace(" ","")
                seq_list.append(new)
    with open(newfile,"w") as g:
        for line in seq_list:
            g.write(line + "\n")
