#!/bin/bash
#SBATCH --mem 10G
#SBATCH --cpus-per-task=5
#SBATCH -o /EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/chemoreceptors/ORs/gff/protein_filtering/logs/tmhmm_%A_%a.out
#SBATCH -C zivnfs
#SBATCH --chdir=/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/chemoreceptors/ORs/gff/protein_filtering

set -e

hostname

FILENAME=$(head -n $SLURM_ARRAY_TASK_ID filenames_tmhmm | tail -1)

mkdir -p biolib_temp_$SLURM_ARRAY_TASK_ID

cd biolib_temp_$SLURM_ARRAY_TASK_ID

cp ../$FILENAME .

biolib run DTU/DeepTMHMM --fasta $FILENAME

cp biolib_results/predicted_topologies.3line ../${FILENAME}.tm
