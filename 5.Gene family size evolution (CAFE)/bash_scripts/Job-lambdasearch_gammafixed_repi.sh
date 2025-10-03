#!/bin/bash

#SBATCH --mail-user=ca3321fu@zedat.fu-berlin.de
#SBATCH --job-name=gfixed                         
#SBATCH --mail-type=none
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=300                                      
#SBATCH --time=50:00:30                                           
#SBATCH --qos=standard
#SBATCH --partition=main,scavenger
#SBATCH --output=slurmMk2/slurm-%j.out

module purge
module load GCC/12.3.0


input_list="M_rep_list.txt"

readarray list < $input_list
Mod2=${list[${SLURM_ARRAY_TASK_ID}-1]}
Mod=$(basename $Mod2)
pre_mod=$(echo $Mod | cut -d "_" -f2 | cut -d "." -f1)
rep=$(dirname $Mod2)

mkdir Runs_for_convergence/${pre_mod}
mkdir Runs_for_convergence/${pre_mod}/k2

/home/ca3321fu/CAFE5/bin/cafe5 -i gene_hog_small_count_eusoc.txt -t tcal_blattodea_tree.tre -y trees_for_optimum/${Mod} -k 1 -p -o Runs_for_convergence/${pre_mod}/k2/${pre_mod}_small_run_k1_rep${rep}