################################################################################

# Author: Cedric Aumont
# Year: 2025

################################################################################


#-------------------------------------------------------------------------------
#-------Introduction -----------------------------------------------------------
#-------------------------------------------------------------------------------

"%ni%" <- Negate("%in%")
len = length
library(wTO)
library(readxl)
library(dplyr)


# set to source file location

#-------------------------------------------------------------------------------
#------- DEG and OR subset -----------------------------------------------------
#-------------------------------------------------------------------------------


spe_list = c("Apac","Cges","Hsjo","Kfla", "Mdar","Mnat","Ncas","PRsim","Rfla")

species_list = c(               # Attention: the order matters: need to be the same as spe_list
  "Anoplotermes_pacificus",
  "Coptotermes_gestroi",
  "Hodotermopsis_sjostedti",
  "Kalotermes_flavicollis",
  "Mastotermes_darwiniensis",
  "Macrotermes_natalensis",
  "Neotermes_castaneus",
  "Prorhinotermes_simplex",
  "Reticulitermes_flavipes")

#------- list of ORs -----------------------------------------------------------

# here we extract the list of ORs

OG112 = read.table("../Data/OG112_spe_gene.txt", col.names = c("spe", "gene"))
OG2gene = read.table("../Data/gene2OG_map.txt",
                     header = T,
                     sep = "\t",
                     dec = ",", 
                     stringsAsFactors = T)
summary(OG2gene)
# formatting species column
spelist = c()
for (line in OG2gene$gene){
  temp = unlist(strsplit(line,"0"))[1]
  if (temp != "Ofor"){
    temp = unlist(strsplit(temp,"O"))[1]}
  temp = unlist(strsplit(temp,"4"))[1]
  spelist = append(spelist,temp)
}
spelist = gsub("Csp", "Csp4", spelist)
spelist = gsub("Bger", "BGER", spelist) # BGER -> Bger of Bger to BGER?
spelist = gsub("KAJ", "Pame", spelist)
OG2gene$spe = spelist
OG2gene$spe = as.factor(OG2gene$spe)
summary(OG2gene)

# formatting node column
nodlist = c()
for (line in OG2gene$gene){
  temp = unlist(strsplit(line,"-"))[1]
  temp = unlist(strsplit(temp,"\\."))[1]
  nodlist = append(nodlist,temp)
}
OG2gene$Node = nodlist
OG2gene$gene = nodlist

OG112 = subset(OG2gene, subset = orthogroup == "OG112", select = c("spe", "gene"))
allOG = subset(OG2gene, select = c("spe", "gene"))


OGclus2 = read.table("../Data/Cluster2_genenames.txt")
colnames(OGclus2) = "gene"
spelist = c()
for (line in OGclus2$gene){
  temp = unlist(strsplit(line,"0"))[1]
  if (temp != "Ofor"){
    temp = unlist(strsplit(temp,"O"))[1]}
  temp = unlist(strsplit(temp,"4"))[1]
  spelist = append(spelist,temp)
}
OGclus2$spe = spelist
OGclus2$spe = as.factor(OGclus2$spe)
summary(OGclus2)


write.table(allOG, "../Result/AllORs_spe_gene.txt",quote = F,sep = "\t",row.names = F)
write.table(OGclus2, "../Result/OGclus2_spe_gene.txt",quote = F,sep = "\t",row.names = F)

#------- DEG all -----------------------------------------------------------
# here we extract the list of DEG

alldeg_wosore = read.table("../TABLES/DEGwosore_genelist_termites.txt", 
                           col.names = c("spe", "gene"))

#------- DEG shared ------------------------------------------------------
# here we want the list of DEG that are shared between 
sHOG_DEGpWORE = read_xlsx(path = "../10.SET/1.DEG/Pairwise_WORE/Result/Results_pDEG_WORE.xlsx",
                          sheet = "sHOG_DEGpWORE-nostat")

# True workers
TW = intersect(subset(sHOG_DEGpWORE, type %in% c("WoBif_5p"))$gene, alldeg_wosore$gene)

# Reproductive (bifurcated)
TR = intersect(subset(sHOG_DEGpWORE, type %in% c("ReBif_5p"))$gene, alldeg_wosore$gene)

# False workers 
FW = intersect(subset(sHOG_DEGpWORE, type %in% c("WoLin_4p"))$gene, alldeg_wosore$gene)

# Reproductive (linear)
FR = intersect(subset(sHOG_DEGpWORE, type %in% c("ReLin_4p"))$gene, alldeg_wosore$gene)




#------- Wo specific ORs -------------------------------------------------------

ORs112_g = data.frame(spe = character(),
                      gene = character(),
                      caste = character())
allOG_g = data.frame(spe = character(),
                     gene = character(),
                     caste = character())
clus2_g = data.frame(spe = character(),
                     gene = character(),
                     caste = character())



for (spe in spe_list){
  species = species_list[match(spe,spe_list)]
  data_WoReGCN = read.csv2(paste0("../Data/data/PairwiseWoReAfterBC_Node_groups_normalised_",species,".csv"))
  
  Wo_ = unique(data_WoReGCN$Group_pval_Phi_tilde[grep(pattern = paste0("g.Network_",species,"_Wo"), data_WoReGCN$Group_pval_Phi_tilde)])
  
  OR_speW = subset(data_WoReGCN,
                   subset = Group_pval_Phi_tilde %in% Wo_ & Node %in% OG112$gene)$Node
  ORall_speW = subset(data_WoReGCN,
                      subset = Group_pval_Phi_tilde %in% Wo_ & Node %in% allOG$gene)$Node
  ORclus2_speW = subset(data_WoReGCN,
                        subset = Group_pval_Phi_tilde %in% Wo_ & Node %in% OGclus2$gene)$Node
  
  
  ORs112_g = bind_rows(ORs112_g, data.frame(spe =rep(spe, len(OR_speW)),
                                            gene = OR_speW,
                                            caste =rep("Wo", len(OR_speW))))
  allOG_g = bind_rows(allOG_g, data.frame(spe =rep(spe, len(ORall_speW)),
                                          gene = ORall_speW,
                                          caste =rep("Wo", len(ORall_speW))))
  clus2_g = bind_rows(clus2_g, data.frame(spe =rep(spe, len(ORclus2_speW)),
                                          gene = ORclus2_speW,
                                          caste =rep("Wo", len(ORclus2_speW))))
  
  
  #------- Re specific ORs -------------------------------------------------------
  
  Re_ = unique(data_WoReGCN$Group_pval_Phi_tilde[grep(pattern = "x_Pr|s_Pr|i_Pr|_Al|_Er|_Nm", data_WoReGCN$Group_pval_Phi_tilde)])
  
  if (spe == "Ncas"){
    Re_ = c("g.Network_Neotermes_castaneus_Al",
            "g.Network_Neotermes_castaneus_Al.Network_Neotermes_castaneus_Er",
            "g.Network_Neotermes_castaneus_Er")
  }
  if (spe == "Kfla"){
    Re_ = c("g.Network_Kalotermes_flavicollis_Al",
            "g.Network_Kalotermes_flavicollis_Al.Network_Kalotermes_flavicollis_Er",
            "g.Network_Kalotermes_flavicollis_Er")
  }
  
  OR_speR = subset(data_WoReGCN,
                   subset = Group_pval_Phi_tilde %in% Re_ & Node %in% OG112$gene)$Node
  ORs112_g = bind_rows(ORs112_g, data.frame(spe =rep(spe, len(OR_speR)),
                                            gene = OR_speR,
                                            caste =rep("Re", len(OR_speR))))
  
  ORall_speR = subset(data_WoReGCN,
                      subset = Group_pval_Phi_tilde %in% Re_ & Node %in% allOG$gene)$Node
  allOG_g = bind_rows(allOG_g, data.frame(spe =rep(spe, len(ORall_speR)),
                                          gene = ORall_speR,
                                          caste =rep("Re", len(ORall_speR))))
  
  ORclus2_speR = subset(data_WoReGCN,
                        subset = Group_pval_Phi_tilde %in% Re_ & Node %in% OGclus2$gene)$Node
  clus2_g = bind_rows(clus2_g, data.frame(spe =rep(spe, len(ORclus2_speR)),
                                          gene = ORclus2_speR,
                                          caste =rep("Re", len(ORclus2_speR))))
}

# write results

write.table(ORs112_g,"../Result/ORs112g_spe_gene_caste_correct.txt",row.names = F,quote = F, sep="\t")
write.table(allOG_g,"../Result/allOGg_spe_gene_caste_correct.txt",row.names = F,quote = F, sep="\t")
write.table(clus2_g,"../Result/clus2g_spe_gene_caste_correct.txt",row.names = F,quote = F, sep="\t")


write.table(TW,"DEG_TW.txt",row.names = F,quote = F, sep="\t")
write.table(FW,"DEG_FW.txt",row.names = F,quote = F, sep="\t")
write.table(TR,"DEG_TR.txt",row.names = F,quote = F, sep="\t")
write.table(FR,"DEG_FR.txt",row.names = F,quote = F, sep="\t")


