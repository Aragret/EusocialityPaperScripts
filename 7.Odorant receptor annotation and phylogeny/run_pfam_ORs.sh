#!/bin/bash
#SBATCH --mem 10G
#SBATCH --cpus-per-task=5
#SBATCH -o /EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/chemoreceptors/ORs/gff/protein_filtering/logs/pfam_%A_%a.out
#SBATCH -C zivnfs
#SBATCH --chdir=/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/chemoreceptors/ORs/gff/protein_filtering

set -e

source /usr/share/modules/init/bash
module use /global/projects/programs/modules/
module load domain/pfam-1.6

FILENAME=$(head -n $SLURM_ARRAY_TASK_ID filenames | tail -1)

pfam_scan.pl -cpu 5 -outfile ${FILENAME}.dom -dir /global/databases/pfam/v37.1/ -fasta $FILENAME

