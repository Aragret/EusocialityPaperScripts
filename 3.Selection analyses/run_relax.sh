#!/bin/bash
#SBATCH --mem 1000M
#SBATCH --cpus-per-task=10
#SBATCH -o logs/relax-%A_%a.out
#SBATCH -C zivnfs
#SBATCH --chdir=/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/relax/without_roaches

##Activating the hyphy module
set -e
source /etc/profile.d/modules.sh 
module unload selection/hyphy
module load selection/hyphy/2.5.41_NOAVX

hostname

#echo $SLURM_ARRAY_TASK_ID 


# all alignment filenames should be in the alignment_files file

maskedfastaFile=$(head -n $SLURM_ARRAY_TASK_ID alignment_files | tail -1)
orthName=${maskedfastaFile%_cds.aln-gb.fas}

#echo ${maskedfastaFile}
#echo ${orthName}

#setting up variables
alignmentFile=/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/cds/og_cds/renamed/${maskedfastaFile}
treeFile=SpeciesTree.txt
outFile=${orthName}-relax.Output

#running hyphy relax
hyphy CPU=10 relax --alignment ${alignmentFile} --tree ${treeFile} --test Foreground --reference Background \
--output ${outFile} > ${outFile}.markdown

