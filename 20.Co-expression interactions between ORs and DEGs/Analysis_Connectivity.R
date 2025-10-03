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


#set wd to source file location
#setwd("~/WORK/Dino_McMahon/Work/Lab work and projects/6. Gene network/OR_to_DEG")


# upload data


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


ORs112_g = read.table("../Data/ORs112g_spe_gene_caste_correct.txt", sep="\t",header = T)
allOG_g = read.table("../Data/allOGg_spe_gene_caste_correct.txt",header = T, sep="\t")
clus2_g = read.table("../Data/clus2g_spe_gene_caste_correct.txt",header = T, sep="\t")

alldeg_wosore = read.table("../Data/DEGwosore_genelist_termites.txt", 
                           col.names = c("spe", "gene"))

sHOG_DEGpWORE = read_xlsx(path = "../Data/Results_pDEG_WORE.xlsx",
                          sheet = "sHOG_DEGpWORE-nostat")

TW = read.table("../Data/DEG_TW.txt", sep="\t",header = T)$x
FW = read.table("../Data/DEG_FW.txt", sep="\t",header = T)$x
TR = read.table("../Data/DEG_TR.txt", sep="\t",header = T)$x
FR = read.table("../Data/DEG_FR.txt", sep="\t",header = T)$x


#-------------------------------------------------------------------------------
#------- Main -----------------------------------------------------
#-------------------------------------------------------------------------------


OR112orall = "all"
# OR112orall = "clus2"
# OR112orall = "112"

df_2 = data.frame(casteNtw = character(),
                  casteOr =  character(),
                  casteDEG =  character(),
                  spe =  character(),
                  res = numeric()
)



for (spe in spe_list){  
  print(spe)
  for (casteNtw in c("Wo", "Re")){
    for (casteOr in c("Wo", "Re")){
      for (casteDEG in c("Wo", "Re")){
        # choose title network
        title = paste0(casteNtw,"NtW_",casteOr,"OR_",casteDEG,"DEG_",spe, "_",OR112orall)
        species = species_list[match(spe,spe_list)]
        if(OR112orall == "all"){
          network = paste0("../Networks/oralldeg_", species, "_",casteNtw,"_nn1.RData")
          load(network)
          if (casteOr == "Wo"){
            ORs = allOG_g$gene[allOG_g$spe == spe & allOG_g$caste == "Wo"]
          } else {
            ORs = allOG_g$gene[allOG_g$spe == spe & allOG_g$caste == "Re"]}
          
        }
        if (OR112orall == "112"){
          network = paste0("../Networks/oralldeg_", species, "_",casteNtw,"_nn1.RData")
          load(network)
          if (casteOr == "Wo"){
            ORs = ORs112_g$gene[ORs112_g$spe == spe & ORs112_g$caste == "Wo"]
          } else {
            ORs = ORs112_g$gene[ORs112_g$spe == spe & ORs112_g$caste == "Re"]}
          
        }
        if (OR112orall == "clus2"){
          network = paste0("../Networks/oralldeg_", species, "_",casteNtw,"_nn1.RData")
          load(network)
          if (casteOr == "Wo"){
            ORs = clus2_g$gene[clus2_g$spe == spe & clus2_g$caste == "Wo"]
          } else {
            ORs = clus2_g$gene[clus2_g$spe == spe & clus2_g$caste == "Re"]}
          
        }
        
        # choose the DEG type
        DEG_spe = sHOG_DEGpWORE[sHOG_DEGpWORE$spe == spe,]
        if (spe %in% c("Apac","Cges","Mdar","Mnat","Rfla")){
          if (casteDEG == "Wo"){onto_deg = TW
          } else {onto_deg = TR}
        } 
        if (spe %in% c("Hsjo","Kfla","Ncas","PRsim")){
          if (casteDEG == "Wo"){onto_deg = FW
          } else {onto_deg = FR}
        } 
        DEG_bia = intersect(subset(DEG_spe, gene %in% onto_deg)$gene, alldeg_wosore$gene)
        
        
        Ntw_OD=rbind.data.frame(
          Network[Network$Node.1 %in% ORs & Network$Node.2 %in% DEG_bia,],
          Network[Network$Node.2 %in% ORs & Network$Node.1 %in% DEG_bia,]
        )
        checkntw = Ntw_OD[abs(Ntw_OD$wTO)>.33,] 
        if (dim(checkntw)[1] >2){
          Vis = NetVis(Ntw_OD$Node.1,
                       Ntw_OD$Node.2, pval = NULL,
                       wTO = Ntw_OD$wTO,MakeGroups =F,#walktrap
                       Cluster =F,cutoff = list(kind="Threshold", value = 0.33),
                       smooth.edges = F,
                       shape = list(shape=c("triangle"), names = c(OG112$gene)),
                       layout = "layout_with_gem")
                       #path = paste0("Network_ORtoDEG2/",title ,".html"))
          
          degree_OR = round(sum(Vis$Nodes[Vis$Nodes$id %in% ORs,]$degree)/len(ORs)/len(DEG_bia),2)
          degree_DEG = round(sum(Vis$Nodes[Vis$Nodes$id %in% DEG_bia,]$degree)/len(DEG_bia)/len(ORs),2)
        }else {
          degree_OR = 0
          degree_DEG = 0
        }
        
        df_2 = bind_rows(df_2, data.frame(casteNtw = casteNtw,
                                          casteOr = casteOr,
                                          casteDEG = casteDEG,
                                          spe = spe,
                                          res = degree_OR
        ))
      }}}}

write.table(df_2, paste0("../Result/analysis_2__OG",OR112orall,"_corrected.txt"), quote = F,row.names = F, sep = "\t")
