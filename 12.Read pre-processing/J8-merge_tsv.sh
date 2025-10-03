#!/bin/bash

#SBATCH --mail-user=cedria95@zedat.fu-berlin.de
#SBATCH --job-name=J8-merge
#SBATCH --mail-type=none
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=400
#SBATCH --time=0:30:00
#SBATCH --qos=standard
#SBATCH --partition=main,scavenger


#input_list="../../lists/species_list.txt"
input_list="new_OR_species_list.txt" # used for the count on the gff file counting the new ORs

readarray list < $input_list
R0=${list[${SLURM_ARRAY_TASK_ID}-1]}
SPE=$(basename $R0)


module purge
module load Biopython

mkdir ${SPE}
cp ../2_hisat2-htseq/${SPE}/*count.tsv ${SPE}

echo "Merging tsv files for ${SPE}"

python3 merge_tsv.py -f ${SPE} -o ${SPE}.readcount.tsv

echo "${SPE}.readcount.tsv is created"


