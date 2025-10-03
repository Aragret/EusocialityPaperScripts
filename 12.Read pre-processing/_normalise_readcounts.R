################################################################################
#
#             Title: read count normalisation 
#           Project: Eusociality
#           Authors: C. Aumont
#              Year: 2024
#
################################################################################


# ------------------------------------------------------------------------------
# ----- [0] Description --------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- [0.1] Purpose ----------------------------------------------------------

'
here the raw readcounts are normalised using DESeq2.
This step was not possible on HPC because of some problems with dependencies
'

# ------------------------------------------------------------------------------
# ------ [1] Introduction ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- [1.1] Load packages ----------------------------------------------------

library(dplyr)
library(DESeq2)

# ----- [1.2] Set directory ----------------------------------------------------

setwd("C:/Users/caumont/Documents/WORK/Dino_McMahon/Work/Lab work and projects/6. Gene network/5.Data preparation/Normalisation_step")


# ----- [1.3] Set arguments ----------------------------------------------------
for (SPECIES in c(
  "Anoplotermes_pacificus",
             "Blatta_orientalis",
             "Blattella_germanica",
             "Coptotermes_gestroi",
             "Cryptocercus_meridianus",
             "Cryptocercus_punctulatus",
             "Hodotermopsis_sjostedti",
             "Kalotermes_flavicollis",
             "Macrotermes_natalensis",
             "Mastotermes_darwiniensis",
             "Neotermes_castaneus",
             "Prorhinotermes_simplex",
             "Reticulitermes_flavipes",
             "Zootermopsis_nevadensis"
             )){
  
readcounts_file = paste0("raw_readcounts/",SPECIES,".readcount.tsv")
conditions_file = paste0("../DESEQ2/conditions/",SPECIES,".conditions.txt")

# ----- [1.] Upload read count ------------------------------------------------

readcounts=read.table(readcounts_file,
                      header=TRUE,
                      row.names=1,
                      sep="\t",
                      stringsAsFactors=FALSE)
dim(readcounts)
head(readcounts)

# ----- [1.] Upload conditions -------------------------------------------------

coldata=read.table(conditions_file,
                   header=TRUE,
                   row.names=1, 
                   sep="\t", 
                   stringsAsFactors=FALSE)
dim(coldata)
head(coldata)


readcounts = readcounts %>% relocate(rownames(coldata))
#View(readcounts)

# ------------------------------------------------------------------------------
# ------ [2] DEG ---------------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- [2.1] factors ----------------------------------------------------------

coldata$caste = factor(coldata$caste)
coldata$sex = factor(coldata$sex)

# ----- [2.2] normalisation ------------------------------------------------------------

dds=DESeqDataSetFromMatrix(countData = readcounts,
                           colData = coldata,
                           design = ~ caste) # + sex + caste:sex
dds

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- round(counts(dds, normalized=TRUE),0)
#view(normalized_counts)

write.table(normalized_counts,
            file=paste0( "normalised_readcounts/",SPECIES,".nreadcount.txt"), 
            sep="\t", quote=F, col.names=NA)
}
