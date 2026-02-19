## OR annotation

- Building HMM profile

`mafft --auto ORs_db.fasta > FPDB_db.aln`

`hmmbuild ORs_db.hmm FPDB_db.aln`

- Running BITACORA for each species with `scripts/runBITACORA_copy.sh` (example for apac)
