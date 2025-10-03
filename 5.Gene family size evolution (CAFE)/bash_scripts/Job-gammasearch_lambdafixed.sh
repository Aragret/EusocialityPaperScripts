#!/bin/bash

#SBATCH --mail-user=ca3321fu@zedat.fu-berlin.de
#SBATCH --job-name=k
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=10:40:30
#SBATCH --qos=standard
#SBATCH --partition=main
#SBATCH --output=slurmK373/slurm-%j.out

module purge
module load GCC/12.3.0


input_list="rate_rep.txt"
readarray list < $input_list
Mod2=${list[${SLURM_ARRAY_TASK_ID}-1]}
rate=$(dirname $Mod2) # k1/rep1
rep=$(basename $Mod2) # k1/rep1

mkdir Runs_for_convergence/gamma373
mkdir Runs_for_convergence/gamma373/k${rate}

#/home/ca3321fu/CAFE5/bin/cafe5 -i gene_hog_small_count_eusoc.txt -t tcal_blattodea_tree.tre -p -k ${rate} -o Runs_for_convergence/gamma/k${rate}/gamma_k${rate}_small_run${rep}
/home/ca3321fu/CAFE5/bin/cafe5 -i gene_hog_small_count_eusoc373.txt -t tcal_blattodea_tree.tre -p -k ${rate} -o Runs_for_convergence/gamma373/k${rate}/gamma_k${rate}_small_run${rep}
