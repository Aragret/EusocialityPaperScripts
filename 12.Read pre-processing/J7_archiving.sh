#!/bin/bash

#SBATCH --mail-user=cedria95@zedat.fu-berlin.de                   # replace with your own address
#SBATCH --job-name=archiving                            # replace the name 
#SBATCH --mail-type=none
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=15G                                        # replace with amount suitable for your job
#SBATCH --time=6:00:00                                           # replace with amount suitable for your job
#SBATCH --qos=standard
#SBATCH --partition=main,scavenger


input_list="/scratch/cedria95/gene_network_analyses/lists/Anoplotermes_pacificus_sample_list.txt"

readarray list < $input_list
R0=${list[${SLURM_ARRAY_TASK_ID}-1]}
R=$(basename $R0)


module purge
module add SAMtools/1.17-GCC-12.2.0

#samtools view -b -S ${R}.sorted.sam > ${R}.sorted.bam
samtools view -h ${R}.sorted.bam > ${R}.sorted.sam
