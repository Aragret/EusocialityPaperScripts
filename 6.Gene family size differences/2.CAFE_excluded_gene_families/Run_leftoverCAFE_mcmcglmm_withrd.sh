#!/bin/bash

#SBATCH --mail-user=cedria95@zedat.fu-berlin.de                   # replace with your own address
#SBATCH --job-name=gcMCMC                          # replace the name 
#SBATCH --mail-type=none
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=360                                        # replace with amount suitable for your job
#SBATCH --time=01:05:30                                           # replace with amount suitable for your job
#SBATCH --qos=standard
#SBATCH --partition=main,scavenger

module load R/4.2.2-foss-2022b

dos2unix data/leftover_list.txt
input_list="data/leftover_list.txt" 

readarray list <$input_list
hog0=${list[${SLURM_ARRAY_TASK_ID} - 1]}
hog=$(basename $hog0)


Rscript leftoverCAFE_MCMCglmm_withrd.R $hog

