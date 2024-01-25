#!/bin/bash
#SBATCH --mem 1000M
#SBATCH --cpus-per-task=4
#SBATCH --array=1-N # job array of the size 1 to the number of alignments
#SBATCH -o absrel-%A-%a.out
#SBATCH --chdir=/working/directory

##Activating the hyphy module
set -e
source /etc/profile.d/modules.sh 
module unload selection/hyphy
module load selection/hyphy/2.5.41_NOAVX

alignment_files=/path/to/file/with/alignment/names

# getting each alignment file and truncating the file name to get orthogroup name
input_alignment=$(head -n $SLURM_ARRAY_TASK_ID $alignment_files | tail -1)
orthName=${input_alignment%_cds.aln-gb.fas}

echo ${orthName}

#setting up variables
alignmentFile=/path/to/${input_alignment}
treeFile=/path/to/labelled.tree
outFile=${orthName}-absrel.Output

#running hyphy relax
hyphy CPU=4 absrel --alignment ${alignmentFile} --tree ${treeFile} --test Foreground \
--output ${outFile} > ${outFile}.markdown
