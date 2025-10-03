#!/bin/bash

#SBATCH --mail-user=ca3321fu@zedat.fu-berlin.de
#SBATCH --job-name=k3_l1_large
#SBATCH --mail-type=none
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1500
#SBATCH --time=20:00:30
#SBATCH --qos=standard
#SBATCH --partition=main
#SBATCH --output=slurmK3_L1_large/slurm-%j.out

module purge
module load GCC/12.3.0


rep=${SLURM_ARRAY_TASK_ID}

mkdir Runs_for_convergence/K1_L3

/home/ca3321fu/CAFE5/bin/cafe5 -i gene_hog_large_count_eusoc.txt -t tcal_blattodea_tree.tre -p -k 3 -l 0.24314831349735 -o Runs_for_convergence/K1_L3/gamma_k3_median_large_run${rep}
