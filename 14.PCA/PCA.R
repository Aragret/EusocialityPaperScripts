################################################################################

# Author: Cedric Aumont
# Year: 2023

# Project: PCA for gene network

################################################################################


#--------------------------------- Aim -----------------------------------------

'
The aim is to identify if caste cluster based on the gene expression for the
different species
'


#-------------------------------------------------------------------------------
#---------------------------------- Directory ----------------------------------
#-------------------------------------------------------------------------------

setwd("~/WORK/Dino_McMahon/Work/Lab work and projects/6. Gene network/4.PCA/All_PCAs")


#-------------------------------------------------------------------------------
#---------------------------------- Library ------------------------------------
#-------------------------------------------------------------------------------

library("FactoMineR")
library("factoextra")


#-------------------------------------------------------------------------------
#---------------------------------- Function -----------------------------------
#-------------------------------------------------------------------------------


do_pca = function(Species, trait, dimA=1, dimB=2){
data0 = read.table(paste0("data/", Species, ".nreadcount.txt"),
                  header = T)
gene_names = rownames(data0)
data1=data0
tdata = t(data1)
colnames(tdata) = gene_names
caste = t(data.frame(strsplit(row.names(tdata),split = "_"))[2,])
sex = t(data.frame(strsplit(row.names(tdata),split = "_"))[3,])
rep = t(data.frame(strsplit(row.names(tdata),split = "_"))[4,])

colnames(caste) = "caste"
colnames(sex) = "sex"
colnames(rep) = "rep"
row.names(caste) = row.names(tdata)
row.names(sex) = row.names(tdata)
row.names(rep) = row.names(tdata)
data = cbind.data.frame(tdata,caste,sex,rep)
head(data[c(1:6),c(1:6)])
title_pca = paste0("PCA - ",
                   strsplit(Species, "_")[[1]][1],
                   " ",
                   strsplit(Species, "_")[[1]][2])

res.pca = PCA(data[,-c(ncol(data)-2,ncol(data)-1,ncol(data))],
              graph = FALSE,
              scale.unit = F)
if(trait == "caste"){
  RES = fviz_pca_ind(res.pca,  axes = c(dimA,dimB),
                     label="none", habillage=as.factor(data$caste),
                     addEllipses=TRUE, ellipse.level=0.7, title=title_pca)
}
if(trait == "sex"){
  RES = fviz_pca_ind(res.pca,  axes = c(dimA, dimB),
                     label="none", habillage=as.factor(data$sex),
                     addEllipses=TRUE, ellipse.level=0.7, title=title_pca)
}
if(trait == "rep"){
  RES = fviz_pca_ind(res.pca,  axes = c(dimA, dimB),
                     label="none", habillage=as.factor(data$rep),
                     addEllipses=TRUE, ellipse.level=0.7, title=title_pca)
}
pdf(file= paste0("results/",trait,"_",Species,"_",dimA,"-",dimB,".pdf"))
plot(RES)
dev.off()
}


#-------------------------------------------------------------------------------
#---------------------------------- MAIN -----------------------------------
#-------------------------------------------------------------------------------

trait_list = c("caste" ,"sex", "rep") 
for (trait in trait_list){
do_pca("Anoplotermes_pacificus", trait)
do_pca("Blatta_orientalis",trait, 1,2)
do_pca("Blattella_germanica",trait, 1,2)
do_pca("Cryptocercus_meridianus",trait,1,2)
do_pca("Cryptocercus_punctulatus",trait,1,2)
do_pca("Mastotermes_darwiniensis",trait,1,2)
do_pca("Macrotermes_natalensis",trait,1,2)
do_pca("Neotermes_castaneus",trait,1,2)
do_pca("Kalotermes_flavicollis",trait,1,2)
do_pca("Hodotermopsis_sjostedti",trait,1,2)
do_pca("Zootermopsis_nevadensis",trait,1,2)
do_pca("Coptotermes_gestroi",trait,1,2)
do_pca("Reticulitermes_flavipes",trait,1,2)
do_pca("Prorhinotermes_simplex",trait,1,2)
}
