#!/bin/bash
#SBATCH --mem 1000M
#SBATCH --cpus-per-task=4
#SBATCH -o /EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/absrel/all_branches/absrel-%A-%a.out
#SBATCH -C zivnfs
#SBATCH --chdir=/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/absrel/all_branches

##Activating the hyphy module
set -e
source /etc/profile.d/modules.sh 
module unload selection/hyphy
module load selection/hyphy/2.5.41_NOAVX

hostname

#echo $SLURM_ARRAY_TASK_ID 

# $SLURM_ARRAY_TASK_ID
maskedfastaFile=$(head -n $SLURM_ARRAY_TASK_ID alignment_files | tail -1)
orthName=${maskedfastaFile%_cds.aln-gb.fas}

#echo ${maskedfastaFile}
echo ${orthName}

#setting up variables
alignmentFile=/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/cds/og_cds/renamed/${maskedfastaFile}
treeFile=/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/codeml/SpeciesTree_mod.txt
outFile=${orthName}-absrel.Output

#running hyphy relax
hyphy CPU=4 absrel --alignment ${alignmentFile} --tree ${treeFile} \
--output ${outFile} > ${outFile}.markdown

