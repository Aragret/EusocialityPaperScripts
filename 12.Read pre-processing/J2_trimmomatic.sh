#!/bin/bash

#SBATCH --mail-user=ca3321fu@zedat.fu-berlin.de                   # replace with your own address
#SBATCH --job-name=J2-Trimmomatic                            # replace the name 
#SBATCH --mail-type=end
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=600                                        # replace with amount suitable for your job
#SBATCH --time=3:00:00                                           # replace with amount suitable for your job
#SBATCH --qos=standard

module purge
module add Trimmomatic/0.39-Java-11

input_list="lists/Zootermopsis_nevadensis_sample_list.txt"

readarray list < $input_list
Ra=${list[${SLURM_ARRAY_TASK_ID}-1]}
R=$(basename $Ra)


java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE brain_transcriptomic_data/${R}_R1.fastq.gz brain_transcriptomic_data/${R}_R2.fastq.gz -baseout Trimmed/${R}.fq.gz ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# note that the Trimmed files have been removed because of lake of space on HPC