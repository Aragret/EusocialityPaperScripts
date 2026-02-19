#!/bin/bash
#SBATCH --mem 4G
#SBATCH --cpus-per-task=1
#SBATCH -o prank/logs/prank-%A_%a.out
#SBATCH --chdir=/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/proteomes/OrthoFinder/Results_Dec13_1/Phylogenetic_Hierarchical_Orthogroups/OtherHOGs_dir/
#SBATCH -C zivnfs

set -e

hostname

# loading prank

source /etc/profile.d/modules.sh

module load msa/prank/170427

# to_run is a file containing filenames of all orthogroups that need to be aligned

OG=$(head -n $SLURM_ARRAY_TASK_ID to_run | tail -1)
orthName=${OG%.fa}

echo ${orthName}

outFile=prank/${orthName}-prank.out

prank -d=${OG} -o=${outFile} -F -once

