#!/bin/bash
#SBATCH --mem 500M
#SBATCH --cpus-per-task=1
#SBATCH --array=1-N # job array of the size 1 to the number of files to align
#SBATCH -o prank-%A_%a.out
#SBATCH --chdir=/working/directory

set -e

# loading prank
source /etc/profile.d/modules.sh
module load msa/prank/170427

# get each orthogroup file name
SingleCopyHOGs = /path/to/file/with/OGfilenames
OG=$(head -n $SLURM_ARRAY_TASK_ID $SingleCopyHOGs | tail -1)

# extract the orthogroup ID
orthName=${OG%.fa}

echo ${orthName}

# setting up variables
outFile=${orthName}-prank.out

prank -d=${OG} -o=${outFile} -F -once
