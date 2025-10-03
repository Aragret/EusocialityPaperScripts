#!/bin/bash

#SBATCH --mail-user=ca3321fu@zedat.fu-berlin.de                   # replace with your own address
#SBATCH --job-name=hisat-built                            # replace the name 
#SBATCH --mail-type=end
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5g                                        # replace with amount suitable for your job
#SBATCH --time=1:30:00                                           # replace with amount suitable for your job
#SBATCH --qos=standard

module purge
module load HISAT2/2.2.1-gompi-2022a

#input_list="lists/species_list.txt"

#readarray list < $input_list
#SPE0=${list[${SLURM_ARRAY_TASK_ID}-1]}
#SPE=$(basename $SPE0)

#mkdir data_preparation/HISAT2_Cufflinks/${SPE}
#hisat2-build -f /scratch/ca3321fu/Termite_Genomes/all_genomes/${SPE}_genome.fna ${SPE}
#mv ${SPE}.* data_preparation/HISAT2_Cufflinks/${SPE}

# for Blattella germanica
#mkdir data_preparation/HISAT2_Cufflinks/Blattella_germanica
#hisat2-build -f /scratch/ca3321fu/Termite_Genomes/Other_genomes/Blattella_germanica_genome.fna Blattella_germanica
#mv Blattella_germanica.* data_preparation/HISAT2_Cufflinks/Blattella_germanica


# for Cryptocercus punctulatus
mkdir data_preparation/HISAT2_Cufflinks/Cryptocercus_punctulatus
hisat2-build -f /scratch/ca3321fu/Termite_Genomes/Other_genomes/Cryptocercus_punctulatus_genome.fna Cryptocercus_punctulatus
mv Cryptocercus_punctulatus.* data_preparation/HISAT2_Cufflinks/Cryptocercus_punctulatus