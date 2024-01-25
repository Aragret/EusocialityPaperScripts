### Creating alignments

- run_prank.sh

### Trim alignments

- run_gblocks.sh

### Concatenate alignments (for phylogeny)

- The program for concatenation that creates a partitioning scheme: https://github.com/PatrickKueck/FASconCAT-G

- Need to run FASconCAT-G_v1.05.1.pl in the directory where the trimmed alignments are.

- Run with `perl ../FASconCAT-G_v1.05.1.pl -s -l` (might modify output options, take a look at the manual).

- Output files that are needed for IQ-TREE2: FcC_supermatrix.fas (concatenated alignment) and FcC_supermatrix_partition.txt (partitioning scheme, I used `-m MFP+MERGE -p FcC_supermatrix_partition.txt` to take partitions into account by iq-tree)
