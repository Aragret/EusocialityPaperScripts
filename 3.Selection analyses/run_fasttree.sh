#!/bin/bash
#SBATCH --mem 500M
#SBATCH --cpus-per-task=1
#SBATCH -o logs/fasttree-%A_%a.out
#SBATCH -C zivnfs

module load phylo/FastTree/2.1.11

set -e

FILENAME=$(head -n $SLURM_ARRAY_TASK_ID ../to_run | tail -1)

echo $FILENAME

FastTreeMP ../${FILENAME} > ${FILENAME}.tree

