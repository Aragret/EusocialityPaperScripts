## Collecting gff files from bitacora output and creating protein and cds files

The following was done with the intermediate Output files of bitacora

### Copying the gff files from the intermediate output

- gff files for curated gene models found in the proteomes:

`while read i; do cp ${i}/ORs/Intermediate_files/ORs_annot_genes_trimmed.gff3 gff/${i}_ORs_annot_genes_trimmed.gff3; done < species_list`

- gff files for curated gene models found in the genomes:

`while read i; do cp ${i}/ORs/Intermediate_files/ORs_genomic_genes_trimmed.gff3 gff/${i}_ORs_genomic_genes_trimmed.gff3; done < species_list`

- concatenating both files:

`while read i; do cat gff/${i}*.gff3 > gff/${i}_ORs_annot_genomic_genes_trimmed.gff3; done < species_list`

- change the species names in gff to match the genome files:

`cd gff/`

`cp apac_ORs_annot_genomic_genes_trimmed.gff3 Apac_ORs_annot_genomic_genes_trimmed.gff3`

### Modifying the attributes in gff

- Renaming split genes

Script: modify_split_ORsg.py

`while read i; do ./modify_split_ORsg.py ${i}_ORs_annot_genomic_genes_trimmed.gff3 ${i}_ORs_annot_genomic_genes_trimmed_renamed.gff3; done < species_list`

- Read mapping program complains about attributes that dont have =

`sed -i 's/blastphmmer/method=blastphmmer/g' *_ORs_annot_genomic_genes_trimmed.gff3`

`sed -i 's/;hmmer;/;method=hmmer;/g' *_renamed.gff3`

`sed -i 's/;-;/;/g' *_renamed.gff3`

`sed -i 's/annot/ann=annot/g' *_renamed.gff3`

- Check for bad attributes:

`awk -F'\t' '                                      
NF == 9 {
    n = split($9, a, ";")
    for (i = 1; i <= n; i++) {
        if (a[i] != "" && a[i] !~ /=/) {
            print "Line " NR ": Bad attribute ->", a[i], "\nFull line: " $0
        }
    }
}' *_renamed.gff3`

- Change the names of new genes by adding the species code

`while read i; do sed "s/=ORsg/=${i}ORsg/g" ${i}_ORs_annot_genomic_genes_trimmed_renamed.gff3 > final/${i}_ORs_annot_genomic_genes_trimmed_renamed.gff3; done < species_list`


### Extracting proteins and applying domain filters

- Make species_list file with removing the part `_ORs_annot_genomic_genes_trimmed.gff3` from the file names

- Extract protein sequences:

`while read i; do gffread final/${i}_ORs_annot_genomic_genes_trimmed_renamed.gff3 -g /EBB-ZIVNFS/amikhail/termite_genomes/annotated_assemblies/gff3/${i}_maskedGenome.fna -y protein_filtering/${i}_ORs_annot_genomic_genes_trimmed_renamed.gff3.pep; done < species_list`

- Run pfam scan to annotate protein domains. Script: run_pfam_ORs.sh

- Extract gene IDs with 7tm_6 domain

`cd protein_filtering/`

`while read i; do grep '7tm_6' ${i}_ORs_annot_genomic_genes_trimmed_renamed.gff3.pep.dom | awk '{print $1}' > ${i}_ORdom_genenames; done < ../species_list`

- Run DeepTMHMM

`sbatch --array=1-29%1 run_tmhmm.sh`

`while read i; do grep " | TM" ${i}_ORs_annot_genomic_genes_trimmed_renamed_ORdom.pep.tm | sed 's/ | TM//g' | sed 's/>//g' > ${i}_ORdom_TM_genenames; done < ../species_list`

### Filtering gff

- Take gene names from the previous step and keep only these genes in gff

`while read i; do grep -B 1 -f ../protein_filtering/${i}_ORdom_TM_genenames ${i}_ORs_annot_genomic_genes_trimmed_renamed.gff3 > filtered/${i}_ORs_annot_genomic_genes_trimmed_renamed_filtered.gff3; done < ../species_list`

- Remove isoform information and new genes to remove ORs from the original annotation

`while read i; do sed 's/-R.*//g' ${i}_ORdom_TM_genenames | sed 's/s1.*//g' | sed 's/s2.*//g' | grep -v 'ORsg' > ${i}_ORdom_TM_genenames_noIso; done < ../species_list `

- Keep the genes not matching OR gene IDs from the original annotation file

`while read i; do grep -v -f ${i}_ORdom_TM_genenames_noIso /EBB-ZIVNFS/amikhail/termite_genomes/annotated_assemblies/gff3/${i}_genes.gff3 > ../final/filtered/${i}_genes_OR_removed.gff3; done < ../species_list`

- Merge OR gff with the rest of the proteome

`while read i; do cat ${i}_genes_OR_removed.gff3 ${i}_ORs_annot_genomic_genes_trimmed_renamed_filtered.gff3 > ${i}_genes_ORs.gff3; done < ../../species_list`

- Make protein and cds files for the filtered OR gff

`while read i; do gffread ${i}_ORs_annot_genomic_genes_trimmed_renamed_filtered.gff3 -g /EBB-ZIVNFS/amikhail/termite_genomes/annotated_assemblies/gff3/${i}_maskedGenome.fna -y ${i}_ORs_annot_genomic_genes_trimmed_renamed_filtered.gff3.pep; done < ../../species_list`

`while read i; do gffread ${i}_ORs_annot_genomic_genes_trimmed_renamed_filtered.gff3 -g /EBB-ZIVNFS/amikhail/termite_genomes/annotated_assemblies/gff3/${i}_maskedGenome.fna -w ${i}_ORs_annot_genomic_genes_trimmed_renamed_filtered.gff3.fna; done < ../../species_list`
