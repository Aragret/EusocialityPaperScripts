#!/bin/bash

#SBATCH --mail-user=ca3321fu@zedat.fu-berlin.de                   
#SBATCH --job-name=rename                            
#SBATCH --mail-type=end
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000                                       
#SBATCH --time=0:10:00                                           
#SBATCH --qos=standard

cd /scratch/ca3321fu/brain_transcriptomic_data

input_listA="Raw_names.txt"
input_listB="Working_names.txt"

readarray listA < $input_listA
nameA=${listA[${SLURM_ARRAY_TASK_ID}-1]}
name_a=$(basename $nameA)

readarray listB < $input_listB
nameB=${listB[${SLURM_ARRAY_TASK_ID}-1]}
name_b=$(basename $nameB)

mv $name_a $name_b
