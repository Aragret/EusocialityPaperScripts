#!/bin/bash

#SBATCH --mail-user=cedria95@zedat.fu-berlin.de                   # replace with your own address
#SBATCH --job-name=htseq_Ap                            # replace the name 
#SBATCH --mail-type=none
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=240                                        # replace with amount suitable for your job
#SBATCH --time=8:00:00                                           # replace with amount suitable for your job
#SBATCH --qos=standard
#SBATCH --partition=main,scavenger


input_list="../../../lists/Anoplotermes_pacificus_sample_list.txt"

readarray list < $input_list
R0=${list[${SLURM_ARRAY_TASK_ID}-1]}
R=$(basename $R0)


module purge
module add SAMtools/1.17-GCC-12.2.0

samtools sort -u -n -O sam -T ${R}.sorted -o ${R}.sorted.sam ${R}.sam



module purge
module add Biopython/1.79-foss-2022a

htseq-count -r name -t gene -i ID -s no -c ${R}.count.tsv ${R}.sorted.sam ../../../gff_afterBC/Apac_genes_ORs.gff3
