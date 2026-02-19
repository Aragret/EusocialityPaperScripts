#!/bin/bash
#SBATCH --job-name=ORs_tree
#SBATCH -o mafft.out
#SBATCH -c 10
#SBATCH --mem=20gb
#SBATCH -C zivnfs

module load msa/mafft/7.505
module load msa/trimal/1.4.1

input_seqs="all_ORs.fa"

# mafft --auto --thread 20 $input_seqs > all_ORs.msa

einsi --thread 10 $input_seqs > all_ORs_einsi_outgroup.msa

trimal -automated1 -in all_ORs_einsi_outgroup.msa -out all_ORs_einsi_outgroup_trimal.msa
