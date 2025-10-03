## Read pre-processing

RNA-seq samples were first renamed

J1_rename.sh

The reads were trimmed using Trimmomatic

J2_trimmomatic.sh

The hisat index was build for each species

J4_hisat2-build.sh

The 3 next steps were ran with a specific script for each species.
Here is shown the script for Apac. The script for the other species 
are the same with the difference that the names, abbreviation, and input 
folder of the species was changed.

J5_Ap_hisat2.sh

J6_Ap_htseq.sh

J7_archiving.sh

J5 maps the genes for each samples using the genome index from the 
previous step and create a sam file for each sample.
J6 first sort the sam file and then use htseq-count to produce the 
raw gene count for each sample.
J7 archives the sam files into bam files.

Further the gene_counts are merged into a single readcount table.

J8-merge_tsv.sh and merge_tsv.py

Finally, the rawreadcount tables are normalised using DESeq2.

_normalise_readcounts.R
