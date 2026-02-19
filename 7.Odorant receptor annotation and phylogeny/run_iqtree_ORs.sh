#!/bin/bash
#SBATCH --mem 40G
#SBATCH --cpus-per-task=20
#SBATCH -o iqtree.out
#SBATCH -C zivnfs

module load phylo/iqtree/2.1.3

set -e

alignment=all_ORs_einsi_outgroup_trimal.msa

iqtree2 -s $alignment -m TEST -ntmax 20
