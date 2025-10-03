#!/bin/bash

#SBATCH --mail-user=
#SBATCH --job-name=l_k_r                         
#SBATCH --mail-type=none
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=400                                        
#SBATCH --time=7:00:30                                           
#SBATCH --qos=standard
#SBATCH --partition=main,scavenger
#SBATCH --output=slurmO/slurm-%j.out

module purge
module load GCC/12.3.0

input_list="l_k_rep.txt"
readarray list < $input_list
Mod2=${list[${SLURM_ARRAY_TASK_ID}-1]}
l=$(echo $Mod2 | cut -d "/" -f1) # l1/k1/rep1
k=$(echo $Mod2 | cut -d "/" -f2)
rep=$(echo $Mod2 | cut -d "/" -f3)

mkdir Runs_for_convergence/cafe_optimal/k${k}_l${l}

/home/ca3321fu/CAFE5/bin/cafe5 -i gene_hog_small_count_eusoc.txt -t tcal_blattodea_tree.tre -y optimal_tree_l${l}.tre -k ${k} -p -o Runs_for_convergence/cafe_optimal/k${k}_l${l}/opt_small_run_k${k}_rep${rep}