#!/bin/bash

#SBATCH --mail-user=ca3321fu@zedat.fu-berlin.de                   # replace with your own address
#SBATCH --job-name=hisat2_Ap                            # replace the name 
#SBATCH --mail-type=end
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3000                                        # replace with amount suitable for your job
#SBATCH --time=13:00:00                                           # replace with amount suitable for your job
#SBATCH --qos=standard

module purge
module load HISAT2/2.2.1-gompi-2022a

input_list="/scratch/ca3321fu/gene_network_analyses/lists/Anoplotermes_pacificus_sample_list.txt"

readarray list < $input_list
R0=${list[${SLURM_ARRAY_TASK_ID}-1]}
R=$(basename $R0)

hisat2 Anoplotermes_pacificus -1 /scratch/ca3321fu/gene_network_analyses/Trimmed/${R}_1P.fq.gz -2 /scratch/ca3321fu/gene_network_analyses/Trimmed/${R}_2P.fq.gz -S ${R}.sam

