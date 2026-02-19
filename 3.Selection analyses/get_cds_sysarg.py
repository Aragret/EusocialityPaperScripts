#!/usr/bin/python3
#get CDS sequences based on protein alignents.

import glob
import MultiFastaLib as mfl
import sys
from pathlib import Path

cds_dict = {}
cdsfiles = glob.glob(sys.argv[1]) #species cds fasta files must be created with gffread so that IDs are same as protein IDs in orthofile
for i in cdsfiles:
    species = i.split("/")[-1][:4] # extract species name from file name
    #species = i.split("/")[-1].replace("_cds_wSpecies.fa_wogene.fa","") # extract species name from file name
    cds_dict[species] = mfl.read_fasta(i) # builds dict with species as key and header-seq-pairs as values (?)

protfiles = glob.glob(sys.argv[2])  # single copy ortholog fastas with gene ids 
for j in protfiles:
    prot_dict = {}
    out_dict = {}
    path = Path(sys.argv[3])
    filename = j.split("/")[-1].replace(".fa","_cds.fasta")
    outfile = Path.joinpath(path, filename)
    with open(j, "r") as f:
        for line in f:
            if ">" in line:
                gene = line.strip().split()[0][1:]
                species = gene[:4]
                prot_dict.setdefault(species, []).append(gene)
    for species in prot_dict:
        for genes in prot_dict.get(species):
            out_dict[genes] = cds_dict[species][genes]

    mfl.write_fasta(out_dict, outfile)


