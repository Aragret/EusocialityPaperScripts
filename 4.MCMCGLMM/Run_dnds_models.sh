#!/bin/bash

#SBATCH --mail-user=cedria95@zedat.fu-berlin.de                   # replace with your own address
#SBATCH --job-name=dndsSCO                          # replace the name 
#SBATCH --mail-type=none
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=460                                        # replace with amount suitable for your job
#SBATCH --time=5-05:05:30                                           # replace with amount suitable for your job
#SBATCH --qos=standard
#SBATCH --partition=main

module load R/4.2.2-foss-2022b

k=${SLURM_ARRAY_TASK_ID}

Rscript MCMCglmm_SCOdnds.R $k
