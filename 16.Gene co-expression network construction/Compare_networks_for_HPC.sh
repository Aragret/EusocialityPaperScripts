#!/bin/bash

#SBATCH --mail-user=cedria95@zedat.fu-berlin.de
#SBATCH --job-name=codina
#SBATCH --mail-type=none
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=392g                                   # replace with amount suitable for your job
#SBATCH --time=2-08:30:00                                           # replace with amount suitable for your job
#SBATCH --qos=standard
#SBATCH --partition=scavenger,main

module load  R/4.3.2-gfbf-2023a
module load  netCDF/4.9.2-gompi-2023a
module load  GSL/2.7-GCC-12.3.0

input_list=/scratch/cedria95/gene_network_analyses/lists/species_list.txt 

readarray list <$input_list
sc=${list[${SLURM_ARRAY_TASK_ID} - 1]}

spe=$(basename $sc)




Rscript Compare_networks_for_HPC.R $spe
