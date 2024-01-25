#!/bin/bash
#SBATCH --mem 2G
#SBATCH --cpus-per-task=5 # might just need one, didn't test
#SBATCH --array=1-N # job array of the size 1 to the number of alignments
#SBATCH -o gblocks-%A-%a.out
#SBATCH --chdir=/working/directory

set -e

FILENAME=$(head -n $SLURM_ARRAY_TASK_ID filenames | tail -1)

# -t=c for codons, -t=p for AAs
/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/src/Gblocks_0.91b/Gblocks $FILENAME -t=c -b5=h
