################################################################################

# Author: Cedric Aumont
# Year: 2023

# Project: PCA for gene network across species

################################################################################


#--------------------------------- Aim -----------------------------------------

'
The aim is to identify if caste cluster based on the gene expression with a single PCA
First, we look at how te samples cluster together when no transformation is applied
then we produce subset of data for caste comparison.
The dataset are averaged in turn by caste and then by sex.
'

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------Libraries------------------------------------------------

library(vegan)
library("FactoMineR")
library("factoextra")
library(SuperExactTest)
"%ni%" = Negate("%in%")
library(readxl)
library(dplyr)

#-------------------------------Work directory------------------------------------------------

setwd("C:/Users/caumont/Documents/WORK/Dino_McMahon/Work/Lab work and projects/6. Gene network/4.PCA/Full_PCA/")


#-------------------------------------------------------------------------------
# ------------------------------Function ---------------------------------------
#-------------------------------------------------------------------------------


do_full_pca = function(data,trait, dimA=1, dimB=2){
  data1=data[,-1]
  tdata = t(data1)
  colnames(tdata) = data$hog_count
  species = t(data.frame(strsplit(row.names(tdata),split = "_"))[1,])
  caste = t(data.frame(strsplit(row.names(tdata),split = "_"))[2,])
  sex = t(data.frame(strsplit(row.names(tdata),split = "_"))[3,])
  colnames(species) = "species"
  colnames(caste) = "caste"
  colnames(sex) = "sex"
  row.names(species) = row.names(tdata)
  row.names(caste) = row.names(tdata)
  row.names(sex) = row.names(tdata)
  data = cbind.data.frame(tdata,species,caste,sex)
  title_pca=""
  res.pca = PCA(data[,-c(ncol(data)-2,ncol(data)-1,ncol(data))],
                graph = FALSE,
                scale.unit = F)
  if(trait == "species"){
    RES = fviz_pca_ind(res.pca,  axes = c(dimA,dimB),
                       label="none", habillage=as.factor(data$species),
                       addEllipses=TRUE, ellipse.level=0.7, title=title_pca)
  }
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
  RES
}
# -------------- normalised ----------------------------------------------------
do_full_norm_pca = function(data,trait, dimA=1, dimB=2){
  data_scaled = data[1]
  spelist = unique(unlist(strsplit(colnames(data[-1]),"_",fixed = T)))
  spelist = spelist[spelist %ni% c("Wo","Re","Re2","So", "Fe", "Ma")]
  for (sp in spelist){
    data_mean = subset(data,
                       select = c("hog",colnames(data)[grep(sp, colnames(data))]))
    row.names(data_mean) = data_mean$hog
    #M = mean(rowMeans(data_mean[-1]))
    #Ncol = dim(data_mean[-1])[2]
    #center_vector = rep(M,Ncol)
    M = (rowMeans(data_mean[-1]))
    center_vector = (M)
    data_scaled_temp = data.frame(scale(t(data_mean[-1]), center = center_vector, scale = center_vector))
    data_scaled_temp = data.frame(t(data_scaled_temp))
    data_scaled_temp$hog  =row.names(data_scaled_temp)
    data_scaled = left_join(data_scaled,data_scaled_temp, by= "hog")
      
  }
  
  
  
  data1=data_scaled[-1]
  tdata = t(data1)
  colnames(tdata) = data$hog_count
  species = t(data.frame(strsplit(row.names(tdata),split = "_"))[1,])
  caste = t(data.frame(strsplit(row.names(tdata),split = "_"))[2,])
  sex = t(data.frame(strsplit(row.names(tdata),split = "_"))[3,])
  colnames(species) = "species"
  colnames(caste) = "caste"
  colnames(sex) = "sex"
  row.names(species) = row.names(tdata)
  row.names(caste) = row.names(tdata)
  row.names(sex) = row.names(tdata)
  data = cbind.data.frame(tdata,species,caste,sex)
  title_pca=""
  res.pca = PCA(data[,-c(ncol(data)-2,ncol(data)-1,ncol(data))],
                graph = FALSE,
                scale.unit = F)
  if(trait == "species"){
    RES = fviz_pca_ind(res.pca,  axes = c(dimA,dimB),
                       label="none", habillage=as.factor(data$species),
                       addEllipses=TRUE, ellipse.level=0.7, title=title_pca)
  }
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
  RES
}
# -------------- return common hogs -------------------------------------------- 
return_common_sog_hogs = function(res, degree, sog){
  deg_res =res$Elements[res$Degree >= degree]
  hog_list =c()
  for (elem in deg_res){
    hog_line = unlist(strsplit(elem,", "))
    for (hog in hog_line){
      if (hog %ni% hog_list){
        hog_list = append(hog_list, hog)
      }  
    }
  }
  unihog = unique(hog_list)
  return(unihog[unihog %in% sog])
}


#-------------------------------------------------------------------------------
#-------------- Dataset preparation --------------------------------------------
#-------------------------------------------------------------------------------

#------- gene to hog file ------------------------------------------------------

gene2hog = read_xlsx(path = "../../TABLES/Tables_SA.xlsx",
                     sheet = "S2_gene_to_hog",
                     skip = 2)


#------- list of single copy orthologues ---------------------------------------
#
SOG_Comp_13_species = data.frame(
  hog=gsub("N0.","",read.table("data/SOG_Comp_13_species.txt",header = T)$HOG))

SOG_All_9_termites = data.frame(
  hog=gsub("N0.","",read.table("data/SOG_All_9_species.txt",header = T)$HOG))

SOG_Lin_termites = data.frame(
  hog=gsub("N0.","",read.table("data/SOG_Lin_species.txt",header = T)$HOG))

SOG_Bif_termites = data.frame(
  hog=gsub("N0.","",read.table("data/SOG_Bif_species.txt",header = T)$HOG))


# ------- creating working dataset ---------------------------------------------

species_list_complete = c(               
  "Anoplotermes_pacificus",
  "Coptotermes_gestroi",
  "Hodotermopsis_sjostedti",
  "Kalotermes_flavicollis",
  "Mastotermes_darwiniensis",
  "Macrotermes_natalensis",
  "Neotermes_castaneus",
  "Prorhinotermes_simplex",
  "Reticulitermes_flavipes",
  "Cryptocercus_meridianus",
  "Cryptocercus_punctulatus",
  "Blattella_germanica",
  "Blatta_orientalis")



species_list_All = c(
  "Anoplotermes_pacificus",
  "Macrotermes_natalensis",
  "Coptotermes_gestroi",
  "Reticulitermes_flavipes",
  "Prorhinotermes_simplex",
  "Neotermes_castaneus",
  "Kalotermes_flavicollis",
  "Hodotermopsis_sjostedti",
  "Mastotermes_darwiniensis"
)

species_list_Lin = c(
  "Prorhinotermes_simplex",
  "Neotermes_castaneus",
  "Kalotermes_flavicollis",
  "Hodotermopsis_sjostedti"
)
species_list_Bif = c(
  "Anoplotermes_pacificus",
  "Macrotermes_natalensis",
  "Coptotermes_gestroi",
  "Reticulitermes_flavipes",
  "Mastotermes_darwiniensis"
)

set_list = c("All","Lin", "Bif","Comp")
for (data_set in set_list){
  if (data_set == "Comp"){
    fin_df = SOG_Comp_13_species
    species_list = species_list_complete
  }
  if (data_set == "All"){
    fin_df = SOG_All_9_termites
    species_list = species_list_All
  }
  if (data_set == "Lin"){
    fin_df = SOG_Lin_termites
    species_list = species_list_Lin
  }
  if (data_set == "Bif"){
    fin_df = SOG_Bif_termites
    species_list = species_list_Bif
  }
for (species in species_list){
# species = "Anoplotermes_pacificus"
# species = "Neotermes_castaneus"
nreadcount = read.table(paste0("data/normalised_count_afterBC/",
                               species,
                               ".nreadcount.txt"))

nreadcount$gene = rownames(nreadcount)
nreadcount = left_join(nreadcount, gene2hog, by="gene")
if (data_set == "Comp"){
  sogcount = subset(nreadcount, hog %in% SOG_Comp_13_species$hog)
}
if (data_set == "All"){
  sogcount = subset(nreadcount, hog %in% SOG_All_9_termites$hog)
}
if (data_set == "Lin"){
  sogcount = subset(nreadcount, hog %in% SOG_Lin_termites$hog)
}
if (data_set == "Bif"){
  sogcount = subset(nreadcount, hog %in% SOG_Bif_termites$hog)
}


list_col = colnames(sogcount)[grep("Wo|Ju", colnames(sogcount))]
sogcount$Wo = rowMeans(sogcount[,colnames(sogcount) %in% list_col])
#if (species %ni% c("Neotermes_castaneus", "Kalotermes_flavicollis")){
  list_col = colnames(sogcount)[grep("Pr|Al|Er|Nm|Ad", colnames(sogcount))]
  sogcount$Re = rowMeans(sogcount[,colnames(sogcount) %in% list_col])
#}else{
#  list_col = colnames(sogcount)[grep("Al", colnames(sogcount))]
#  sogcount$Re = rowMeans(sogcount[,colnames(sogcount) %in% list_col])
#  list_col = colnames(sogcount)[grep("Er", colnames(sogcount))]
#  sogcount$Re2 = rowMeans(sogcount[,colnames(sogcount) %in% list_col])
#}
list_col = colnames(sogcount)[grep("So", colnames(sogcount))] # added for So
sogcount$So = rowMeans(sogcount[,colnames(sogcount) %in% list_col]) # added for So

sp = unlist(strsplit(colnames(nreadcount)[1],"_"))[1]

# if ("Re2" %in% colnames(sogcount)){
#   inter_df = subset(sogcount, select = c(hog, Wo, Re, Re2))
#   colnames(inter_df) = c("hog", paste0(sp,"_Wo"),paste0(sp,"_Re")
#                          ,paste0(sp,"_Re2"))
# }else{
#   inter_df = subset(sogcount, select = c(hog, Wo, Re))
#   colnames(inter_df) = c("hog", paste0(sp,"_Wo"),paste0(sp,"_Re"))
# }
# changed for So 
if ("Re2" %in% colnames(sogcount)){
  inter_df = subset(sogcount, select = c(hog, Wo, Re, Re2, So))
  colnames(inter_df) = c("hog", paste0(sp,"_Wo"),paste0(sp,"_Re")
                         ,paste0(sp,"_Re2"),paste0(sp,"_So"))
}else{
  inter_df = subset(sogcount, select = c(hog, Wo, Re, So))
  colnames(inter_df) = c("hog", paste0(sp,"_Wo"),paste0(sp,"_Re"),paste0(sp,"_So"))
}

fin_df = left_join(fin_df, inter_df, by="hog")
}
saveRDS(fin_df, paste0("data/data_for_",data_set,"_species_WoSoRe.RDS"))
}


#-------------------------------------------------------------------------------
#-------------- Preparing degree selection--------------------------------------
#-------------------------------------------------------------------------------

# ----------- Wo_Re ------------------------------------------------------------

datas_Re = readRDS("../../10.SET/1.DEG/Pairwise_WORE/Result/data_biased_Re.rds")
datas_Wo = readRDS("../../10.SET/1.DEG/Pairwise_WORE/Result/data_biased_Wo.rds")



# -- 

######
##
# Comp
data_Re_termites = datas_Re
data_Wo_termites = datas_Wo
Result=supertest(data_Re_termites)
res=data.frame(summary(Result)$Table)
hog_Re_13 = return_common_sog_hogs(res = res, degree = 13, sog = SOG_Comp_13_species$hog)

hog_Re_10 = return_common_sog_hogs(res = res, degree = 10, sog = SOG_Comp_13_species$hog)
hog_Re_9 = return_common_sog_hogs(res = res, degree = 9, sog = SOG_Comp_13_species$hog)
hog_Re_8 = return_common_sog_hogs(res = res, degree = 8, sog = SOG_Comp_13_species$hog)
hog_Re_7 = return_common_sog_hogs(res = res, degree = 7, sog = SOG_Comp_13_species$hog)
hog_Re_6 = return_common_sog_hogs(res = res, degree = 6, sog = SOG_Comp_13_species$hog)
hog_Re_5 = return_common_sog_hogs(res = res, degree = 5, sog = SOG_Comp_13_species$hog)
hog_Re_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_Comp_13_species$hog)
hog_Re_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_Comp_13_species$hog)
hog_Re_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_Comp_13_species$hog)
hog_Re_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_Comp_13_species$hog)

Result=supertest(data_Wo_termites)
res=data.frame(summary(Result)$Table)
hog_Wo_13 = return_common_sog_hogs(res = res, degree = 13, sog = SOG_Comp_13_species$hog)
hog_Wo_10 = return_common_sog_hogs(res = res, degree = 10, sog = SOG_Comp_13_species$hog)
hog_Wo_9 = return_common_sog_hogs(res = res, degree = 9, sog = SOG_Comp_13_species$hog)
hog_Wo_8 = return_common_sog_hogs(res = res, degree = 8, sog = SOG_Comp_13_species$hog)
hog_Wo_7 = return_common_sog_hogs(res = res, degree = 7, sog = SOG_Comp_13_species$hog)
hog_Wo_6 = return_common_sog_hogs(res = res, degree = 6, sog = SOG_Comp_13_species$hog)
hog_Wo_5 = return_common_sog_hogs(res = res, degree = 5, sog = SOG_Comp_13_species$hog)
hog_Wo_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_Comp_13_species$hog)
hog_Wo_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_Comp_13_species$hog)
hog_Wo_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_Comp_13_species$hog)
hog_Wo_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_Comp_13_species$hog)

Hog_Com_13p = c(hog_Re_13,hog_Wo_13)
Hog_Com_10p = c(hog_Re_10,hog_Wo_10)
Hog_Com_9p = c(hog_Re_9,hog_Wo_9)
Hog_Com_8p = c(hog_Re_8,hog_Wo_8)
Hog_Com_7p = c(hog_Re_7,hog_Wo_7)
Hog_Com_6p = c(hog_Re_6,hog_Wo_6)
Hog_Com_5p = c(hog_Re_5,hog_Wo_5)
Hog_Com_4p = c(hog_Re_4,hog_Wo_4)
Hog_Com_3p = c(hog_Re_3,hog_Wo_3)
Hog_Com_2p = c(hog_Re_2,hog_Wo_2)
Hog_Com_1p = c(hog_Re_1,hog_Wo_1)

#data_comp = readRDS("data/data_for_Comp_termites_withSo.RDS")
data_comp = readRDS("data/data_for_Comp_species_WoSoRe.RDS")


pdf("norm_PCA_Comp_WoRe_13sp_7p_shared_genes_species.pdf")
do_full_norm_pca(subset(data_comp,subset = hog %in% Hog_Com_7p,
                        select = colnames(data_comp) %ni% c(colnames(data_comp)[grep("So",colnames(data_comp))])),
                 trait = "species")
dev.off()
dim(subset(data_comp,subset = hog %in% Hog_Com_7p,
           select = colnames(data_comp) %ni% c(colnames(data_comp)[grep("So",colnames(data_comp))],
                                               "Kf_Re2","Nc_Re2")))
dimensionSOG = c(3464,2247,1088,480,190,66,17)

do_full_pca(subset(data_comp,subset = hog %in% Hog_Com_1p,
                   select = colnames(data_comp) %ni% c("Ap_So","Cm_So","Bg_So","Bo_So",
                                                       "Cp_So","Kf_Re2","Nc_Re2")),
            trait = "species")

data_comp2 = data_comp
colnames(data_comp2) = c("hog",  "Ap_Wo",  "Ap_Re",  "Ap_So",  "Cg_Wo",  "Cg_Re",  "Cg_So",
                         "Hs_Wo",  "Hs_Re",  "Hs_So", "Kf_Wo",  "Kf_Re",  "Kf_Re2", "Kf_So",
                         "Md_Wo",  "Md_Re",  "Md_So",  "Mn_Wo",  "Mn_Re",  "Mn_So", "Nc_Wo",
                         "Nc_Re",  "Nc_Re2", "Nc_So",  "Ps_Wo" , "Ps_Re" , "Ps_So","Rf_Wo",
                         "Rf_Re"  ,"Rf_So" ,"Cm_Wo" , "Cm_Ad" , "Cm_So" , "Cp_Wo" , "Cp_Ad",
                         "Cp_So", "Bg_Wo" , "Bg_Ad",  "Bg_So","Bo_Wo" ,"Bo_Ad",  "Bo_So") 
do_full_norm_pca(subset(data_comp,subset = hog %in% Hog_Com_1p,
                        select = colnames(data_comp) %ni% c(colnames(data_comp)[grep("So",colnames(data_comp))],
                                                            "Kf_Re2","Nc_Re2","Cm_Wo", "Cm_Re", "Ap_Wo","Ap_Re")),
                 trait = "species")
do_full_norm_pca(subset(data_comp,subset = hog %in% Hog_Com_4p,
                        select = colnames(data_comp) %ni% c(colnames(data_comp)[grep("So",colnames(data_comp))])),
                 trait = "caste")

dim(subset(data_comp,subset = hog %in% Hog_Com_4p,
           select = colnames(data_comp) %ni% c(colnames(data_comp)[grep("So",colnames(data_comp))],
                                               "Kf_Re2","Nc_Re2")))
##
######


data_Re_termites = datas_Re[names(datas_Re) %ni% c("Cmer_Re", "Cpun_Re",
                                                   "BGER_Re", "Bori_Re")]
data_Wo_termites = datas_Wo[names(datas_Wo) %ni% c("Cmer_Wo", "Cpun_Wo", 
                                                   "BGER_Wo", "Bori_Wo")]

Result=supertest(data_Re_termites)
res=data.frame(summary(Result)$Table)
hog_Re_9 = return_common_sog_hogs(res = res, degree = 9, sog = SOG_All_9_termites$hog)
hog_Re_8 = return_common_sog_hogs(res = res, degree = 8, sog = SOG_All_9_termites$hog)
hog_Re_7 = return_common_sog_hogs(res = res, degree = 7, sog = SOG_All_9_termites$hog)
hog_Re_6 = return_common_sog_hogs(res = res, degree = 6, sog = SOG_All_9_termites$hog)
hog_Re_5 = return_common_sog_hogs(res = res, degree = 5, sog = SOG_All_9_termites$hog)
hog_Re_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_All_9_termites$hog)
hog_Re_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_All_9_termites$hog)
hog_Re_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_All_9_termites$hog)
hog_Re_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_All_9_termites$hog)

Result=supertest(data_Wo_termites)
res=data.frame(summary(Result)$Table)
hog_Wo_9 = return_common_sog_hogs(res = res, degree = 9, sog = SOG_All_9_termites$hog)
hog_Wo_8 = return_common_sog_hogs(res = res, degree = 8, sog = SOG_All_9_termites$hog)
hog_Wo_7 = return_common_sog_hogs(res = res, degree = 7, sog = SOG_All_9_termites$hog)
hog_Wo_6 = return_common_sog_hogs(res = res, degree = 6, sog = SOG_All_9_termites$hog)
hog_Wo_5 = return_common_sog_hogs(res = res, degree = 5, sog = SOG_All_9_termites$hog)
hog_Wo_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_All_9_termites$hog)
hog_Wo_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_All_9_termites$hog)
hog_Wo_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_All_9_termites$hog)
hog_Wo_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_All_9_termites$hog)

Hog_All_9p = c(hog_Re_9,hog_Wo_9)
Hog_All_8p = c(hog_Re_8,hog_Wo_8)
Hog_All_7p = c(hog_Re_7,hog_Wo_7)
Hog_All_6p = c(hog_Re_6,hog_Wo_6)
Hog_All_5p = c(hog_Re_5,hog_Wo_5)
Hog_All_4p = c(hog_Re_4,hog_Wo_4)
Hog_All_3p = c(hog_Re_3,hog_Wo_3)
Hog_All_2p = c(hog_Re_2,hog_Wo_2)
Hog_All_1p = c(hog_Re_1,hog_Wo_1)


data_Re_linear = datas_Re[names(datas_Re) %in% c("Hsjo_Re", "Kfla_Re",
                                                 "Ncas_Re", "PRsim_Re")]
data_Wo_linear = datas_Wo[names(datas_Wo) %in% c("Hsjo_Wo", "Kfla_Wo",
                                                 "Ncas_Wo", "PRsim_Wo")]

Result=supertest(data_Re_linear)
res=data.frame(summary(Result)$Table)
hog_Re_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_Lin_termites$hog)
hog_Re_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_Lin_termites$hog)
hog_Re_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_Lin_termites$hog)
hog_Re_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_Lin_termites$hog)

Result=supertest(data_Wo_linear)
res=data.frame(summary(Result)$Table)
hog_Wo_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_Lin_termites$hog)
hog_Wo_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_Lin_termites$hog)
hog_Wo_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_Lin_termites$hog)
hog_Wo_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_Lin_termites$hog)


Hog_Lin_4p = c(hog_Re_4,hog_Wo_4)
Hog_Lin_3p = c(hog_Re_3,hog_Wo_3)
Hog_Lin_2p = c(hog_Re_2,hog_Wo_2)
Hog_Lin_1p = c(hog_Re_1,hog_Wo_1)

data_Wo_bifur = datas_Wo[names(datas_Wo) %in% c("Mdar_Wo", "Rfla_Wo", "Cges_Wo",
                                                "Mnat_Wo", "Apac_Wo")]
data_Re_bifur = datas_Re[names(datas_Re) %in% c("Mdar_Re", "Rfla_Re", "Cges_Re",
                                                "Mnat_Re", "Apac_Re")]

Result=supertest(data_Re_bifur)
res=data.frame(summary(Result)$Table)
hog_Re_5 = return_common_sog_hogs(res = res, degree = 5, sog = SOG_Bif_termites$hog)
hog_Re_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_Bif_termites$hog)
hog_Re_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_Bif_termites$hog)
hog_Re_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_Bif_termites$hog)
hog_Re_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_Bif_termites$hog)

Result=supertest(data_Wo_bifur)
res=data.frame(summary(Result)$Table)
hog_Wo_5 = return_common_sog_hogs(res = res, degree = 5, sog = SOG_Bif_termites$hog)
hog_Wo_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_Bif_termites$hog)
hog_Wo_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_Bif_termites$hog)
hog_Wo_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_Bif_termites$hog)
hog_Wo_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_Bif_termites$hog)

Hog_Bif_5p = c(hog_Re_5,hog_Wo_5)
Hog_Bif_4p = c(hog_Re_4,hog_Wo_4)
Hog_Bif_3p = c(hog_Re_3,hog_Wo_3)
Hog_Bif_2p = c(hog_Re_2,hog_Wo_2)
Hog_Bif_1p = c(hog_Re_1,hog_Wo_1)

library(dplyr)

data_all = readRDS("data/data_for_All_termites.RDS")
data_lin = readRDS("data/data_for_Lin_termites.RDS")
data_bif = readRDS("data/data_for_Bif_termites.RDS")

do_full_pca(subset(data_all,subset = hog %in% Hog_All_1p,
                   select = colnames(data_all) %ni% c("Kf_Re2","Nc_Re2")),
                   trait = "caste")

do_full_norm_pca(subset(data_all,subset = hog %in% Hog_All_1p,
                   select = colnames(data_all) %ni% c("Kf_Re2","Nc_Re2")),
            trait = "caste")
do_full_norm_pca(subset(data_all,subset = hog %in% Hog_All_1p),
                 trait = "species")

do_full_norm_pca(subset(data_all,subset = hog %in% Hog_All_4p),
                 trait = "caste")
#datax = subset(data,subset = hog %in% Hog_All_9p,
 #              select = colnames(data) %ni% c("Kf_Re2","Nc_Re2"))

do_full_pca(subset(data_lin,subset = hog %in% Hog_Lin_1p,
                   select = colnames(data_lin) %ni% c("Kf_Re2","Nc_Re2")),
            trait = "caste")

do_full_norm_pca(subset(data_lin,subset = hog %in% Hog_Lin_1p,
                        select = colnames(data_lin) %ni% c("Kf_Re2","Nc_Re2")),
                 trait = "caste")

do_full_pca(subset(data_bif,subset = hog %in% Hog_Bif_1p,
                   select = colnames(data_bif) %ni% c("Kf_Re2","Nc_Re2")),
            trait = "caste")

do_full_norm_pca(subset(data_bif,subset = hog %in% Hog_Bif_1p,
                        select = colnames(data_bif) %ni% c("Kf_Re2","Nc_Re2")),
                 trait = "caste")

##

# ----------- Wo_So ------------------------------------------------------------

datas_wSo = readRDS("../../10.SET/1.DEG/Pairwise_WOSO/Result/data_biased_So.rds")
datas_sWo = readRDS("../../10.SET/1.DEG/Pairwise_WOSO/Result/data_biased_Wo.rds")


data_wSo_termites = datas_wSo
data_sWo_termites = datas_sWo

Result=supertest(data_wSo_termites)
res=data.frame(summary(Result)$Table)
hog_So_9 = return_common_sog_hogs(res = res, degree = 9, sog = SOG_All_9_termites$hog)
hog_So_8 = return_common_sog_hogs(res = res, degree = 8, sog = SOG_All_9_termites$hog)
hog_So_7 = return_common_sog_hogs(res = res, degree = 7, sog = SOG_All_9_termites$hog)
hog_So_6 = return_common_sog_hogs(res = res, degree = 6, sog = SOG_All_9_termites$hog)
hog_So_5 = return_common_sog_hogs(res = res, degree = 5, sog = SOG_All_9_termites$hog)
hog_So_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_All_9_termites$hog)
hog_So_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_All_9_termites$hog)
hog_So_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_All_9_termites$hog)
hog_So_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_All_9_termites$hog)

Result=supertest(data_sWo_termites)
res=data.frame(summary(Result)$Table)
hog_sWo_9 = return_common_sog_hogs(res = res, degree = 9, sog = SOG_All_9_termites$hog)
hog_sWo_8 = return_common_sog_hogs(res = res, degree = 8, sog = SOG_All_9_termites$hog)
hog_sWo_7 = return_common_sog_hogs(res = res, degree = 7, sog = SOG_All_9_termites$hog)
hog_sWo_6 = return_common_sog_hogs(res = res, degree = 6, sog = SOG_All_9_termites$hog)
hog_sWo_5 = return_common_sog_hogs(res = res, degree = 5, sog = SOG_All_9_termites$hog)
hog_sWo_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_All_9_termites$hog)
hog_sWo_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_All_9_termites$hog)
hog_sWo_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_All_9_termites$hog)
hog_sWo_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_All_9_termites$hog)

Hogsw_All_9p = c(hog_So_9,hog_sWo_9)
Hogsw_All_8p = c(hog_So_8,hog_sWo_8)
Hogsw_All_7p = c(hog_So_7,hog_sWo_7)
Hogsw_All_6p = c(hog_So_6,hog_sWo_6)
Hogsw_All_5p = c(hog_So_5,hog_sWo_5)
Hogsw_All_4p = c(hog_So_4,hog_sWo_4)
Hogsw_All_3p = c(hog_So_3,hog_sWo_3)
Hogsw_All_2p = c(hog_So_2,hog_sWo_2)
Hogsw_All_1p = c(hog_So_1,hog_sWo_1)

#
data_wSo_linear = datas_wSo[names(datas_wSo) %in% c("Hsjo_So", "Kfla_So",
                                                 "Ncas_So", "PRsim_So")]
data_sWo_linear = datas_sWo[names(datas_sWo) %in% c("Hsjo_Wo", "Kfla_Wo",
                                                 "Ncas_Wo", "PRsim_Wo")]

Result=supertest(data_wSo_linear)
res=data.frame(summary(Result)$Table)
hog_So_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_Lin_termites$hog)
hog_So_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_Lin_termites$hog)
hog_So_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_Lin_termites$hog)
hog_So_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_Lin_termites$hog)

Result=supertest(data_sWo_linear)
res=data.frame(summary(Result)$Table)
hog_Wo_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_Lin_termites$hog)
hog_Wo_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_Lin_termites$hog)
hog_Wo_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_Lin_termites$hog)
hog_Wo_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_Lin_termites$hog)


Hogsw_Lin_4p = unique(c(Hog_Lin_4p,hog_Re_4,hog_Wo_4))
Hogsw_Lin_3p = unique(c(Hog_Lin_3p,hog_Re_3,hog_Wo_3))
Hogsw_Lin_2p = unique(c(Hog_Lin_2p,hog_Re_2,hog_Wo_2))
Hogsw_Lin_1p = unique(c(Hog_Lin_1p,hog_Re_1,hog_Wo_1))

data_sWo_bifur = datas_sWo[names(datas_sWo) %in% c("Mdar_Wo", "Rfla_Wo", "Cges_Wo",
                                                "Mnat_Wo")]
data_wSo_bifur = datas_wSo[names(datas_wSo) %in% c("Mdar_So", "Rfla_So", "Cges_So",
                                                "Mnat_So")]

Result=supertest(data_wSo_bifur)
res=data.frame(summary(Result)$Table)
hog_So_5 = return_common_sog_hogs(res = res, degree = 5, sog = SOG_Bif_termites$hog)
hog_So_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_Bif_termites$hog)
hog_So_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_Bif_termites$hog)
hog_So_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_Bif_termites$hog)
hog_So_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_Bif_termites$hog)

Result=supertest(data_sWo_bifur)
res=data.frame(summary(Result)$Table)
hog_Wo_5 = return_common_sog_hogs(res = res, degree = 5, sog = SOG_Bif_termites$hog)
hog_Wo_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_Bif_termites$hog)
hog_Wo_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_Bif_termites$hog)
hog_Wo_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_Bif_termites$hog)
hog_Wo_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_Bif_termites$hog)

Hogsw_Bif_5p = unique(c(Hog_Bif_5p,hog_So_5,hog_Wo_5))
Hogsw_Bif_4p = unique(c(Hog_Bif_4p,hog_So_4,hog_Wo_4))
Hogsw_Bif_3p = unique(c(Hog_Bif_3p,hog_So_3,hog_Wo_3))
Hogsw_Bif_2p = unique(c(Hog_Bif_2p,hog_So_2,hog_Wo_2))
Hogsw_Bif_1p = unique(c(Hog_Bif_1p,hog_So_1,hog_Wo_1))





#data_allso = readRDS("data/data_for_All_termites_withSo.RDS")
data_allso = readRDS("data/data_for_All_species_WoSoRe.RDS")
data_linso = readRDS("data/data_for_Lin_species_WoSoRe.RDS")
data_bifso = readRDS("data/data_for_Bif_species_WoSoRe.RDS")

Hog_All_1p_t = unique(c(Hog_All_1p, Hogsw_All_1p))
Hog_All_2p_t = unique(c(Hog_All_2p, Hogsw_All_2p))
Hog_All_3p_t = unique(c(Hog_All_3p, Hogsw_All_3p))
Hog_All_4p_t = unique(c(Hog_All_4p, Hogsw_All_4p))
Hog_All_5p_t = unique(c(Hog_All_5p, Hogsw_All_5p))
Hog_All_6p_t = unique(c(Hog_All_6p, Hogsw_All_6p))
Hog_All_7p_t = unique(c(Hog_All_7p, Hogsw_All_7p))



do_full_norm_pca(subset(data_allso,subset = hog %in% Hog_All_1p_t,
                   select = colnames(data_allso) %ni% c("Zn_Re","Ap_So","Kf_Re2","Nc_Re2")),
            trait = "caste")

do_full_norm_pca(subset(data_allso,subset = hog %in% Hog_All_6p_t,
                        select = colnames(data_allso) %ni% c("Zn_Re","Ap_So","Kf_Re2","Nc_Re2")),
                 trait = "species")
do_full_norm_pca(subset(data_all,subset = hog %in% Hog_All_4p),
                 trait = "caste")
dim(unique(data_allso))

#




# ----------- So_Re ------------------------------------------------------------

datas_sRe = readRDS("../../10.SET/1.DEG/Pairwise_SORE/Result/data_biased_Re.rds")
datas_rSo = readRDS("../../10.SET/1.DEG/Pairwise_SORE/Result/data_biased_So.rds")


# -- 



data_sRe_termites = datas_sRe[names(datas_sRe) %ni% c("Cmer_Re", "Cpun_Re",
                                                   "BGER_Re", "Bori_Re")]
data_rSo_termites = datas_rSo[names(datas_rSo) %ni% c("Cmer_Wo", "Cpun_Wo", 
                                                   "BGER_Wo", "Bori_Wo")]

Result=supertest(data_sRe_termites)
res=data.frame(summary(Result)$Table)
hog_Re_9 = return_common_sog_hogs(res = res, degree = 9, sog = SOG_All_9_termites$hog)
hog_Re_8 = return_common_sog_hogs(res = res, degree = 8, sog = SOG_All_9_termites$hog)
hog_Re_7 = return_common_sog_hogs(res = res, degree = 7, sog = SOG_All_9_termites$hog)
hog_Re_6 = return_common_sog_hogs(res = res, degree = 6, sog = SOG_All_9_termites$hog)
hog_Re_5 = return_common_sog_hogs(res = res, degree = 5, sog = SOG_All_9_termites$hog)
hog_Re_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_All_9_termites$hog)
hog_Re_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_All_9_termites$hog)
hog_Re_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_All_9_termites$hog)
hog_Re_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_All_9_termites$hog)

Result=supertest(data_rSo_termites)
res=data.frame(summary(Result)$Table)
hog_Wo_9 = return_common_sog_hogs(res = res, degree = 9, sog = SOG_All_9_termites$hog)
hog_Wo_8 = return_common_sog_hogs(res = res, degree = 8, sog = SOG_All_9_termites$hog)
hog_Wo_7 = return_common_sog_hogs(res = res, degree = 7, sog = SOG_All_9_termites$hog)
hog_Wo_6 = return_common_sog_hogs(res = res, degree = 6, sog = SOG_All_9_termites$hog)
hog_Wo_5 = return_common_sog_hogs(res = res, degree = 5, sog = SOG_All_9_termites$hog)
hog_Wo_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_All_9_termites$hog)
hog_Wo_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_All_9_termites$hog)
hog_Wo_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_All_9_termites$hog)
hog_Wo_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_All_9_termites$hog)

Hogsr_All_9p = c(hog_Re_9,hog_Wo_9)
Hogsr_All_8p = c(hog_Re_8,hog_Wo_8)
Hogsr_All_7p = c(hog_Re_7,hog_Wo_7)
Hogsr_All_6p = c(hog_Re_6,hog_Wo_6)
Hogsr_All_5p = c(hog_Re_5,hog_Wo_5)
Hogsr_All_4p = c(hog_Re_4,hog_Wo_4)
Hogsr_All_3p = c(hog_Re_3,hog_Wo_3)
Hogsr_All_2p = c(hog_Re_2,hog_Wo_2)
Hogsr_All_1p = c(hog_Re_1,hog_Wo_1)


data_sRe_linear = datas_sRe[names(datas_sRe) %in% c("Hsjo_Re", "Kfla_Re",
                                                 "Ncas_Re", "PRsim_Re")]
data_rSo_linear = datas_rSo[names(datas_rSo) %in% c("Hsjo_So", "Kfla_So",
                                                 "Ncas_So", "PRsim_So")]

Result=supertest(data_sRe_linear)
res=data.frame(summary(Result)$Table)
hog_Re_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_Lin_termites$hog)
hog_Re_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_Lin_termites$hog)
hog_Re_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_Lin_termites$hog)
hog_Re_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_Lin_termites$hog)

Result=supertest(data_rSo_linear)
res=data.frame(summary(Result)$Table)
hog_So_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_Lin_termites$hog)
hog_So_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_Lin_termites$hog)
hog_So_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_Lin_termites$hog)
hog_So_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_Lin_termites$hog)


Hogsr_Lin_4p = unique(c(Hogsw_Lin_4p,hog_Re_4,hog_So_4))
Hogsr_Lin_3p = unique(c(Hogsw_Lin_3p,hog_Re_3,hog_So_3))
Hogsr_Lin_2p = unique(c(Hogsw_Lin_2p,hog_Re_2,hog_So_2))
Hogsr_Lin_1p = unique(c(Hogsw_Lin_1p,hog_Re_1,hog_So_1))


data_rSo_bifur = datas_rSo[names(datas_rSo) %in% c("Mdar_So", "Rfla_So", "Cges_So",
                                                "Mnat_So")]
data_sRe_bifur = datas_sRe[names(datas_sRe) %in% c("Mdar_Re", "Rfla_Re", "Cges_Re",
                                                "Mnat_Re")]

Result=supertest(data_sRe_bifur)
res=data.frame(summary(Result)$Table)
hog_Re_5 = return_common_sog_hogs(res = res, degree = 5, sog = SOG_Bif_termites$hog)
hog_Re_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_Bif_termites$hog)
hog_Re_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_Bif_termites$hog)
hog_Re_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_Bif_termites$hog)
hog_Re_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_Bif_termites$hog)

Result=supertest(data_rSo_bifur)
res=data.frame(summary(Result)$Table)
hog_So_5 = return_common_sog_hogs(res = res, degree = 5, sog = SOG_Bif_termites$hog)
hog_So_4 = return_common_sog_hogs(res = res, degree = 4, sog = SOG_Bif_termites$hog)
hog_So_3 = return_common_sog_hogs(res = res, degree = 3, sog = SOG_Bif_termites$hog)
hog_So_2 = return_common_sog_hogs(res = res, degree = 2, sog = SOG_Bif_termites$hog)
hog_So_1 = return_common_sog_hogs(res = res, degree = 1, sog = SOG_Bif_termites$hog)

Hog_Bif_5p_tt = unique(c(Hogsw_Bif_5p,hog_Re_5,hog_So_5))
Hog_Bif_4p_tt = unique(c(Hogsw_Bif_4p,hog_Re_4,hog_So_4))
Hog_Bif_3p_tt = unique(c(Hogsw_Bif_3p,hog_Re_3,hog_So_3))
Hog_Bif_2p_tt = unique(c(Hogsw_Bif_2p,hog_Re_2,hog_So_2))
Hog_Bif_1p_tt = unique(c(Hogsw_Bif_1p,hog_Re_1,hog_So_1))


Hog_All_1p_tt = unique(c(Hog_All_1p_t, Hogsr_All_1p))
Hog_All_2p_tt = unique(c(Hog_All_2p_t, Hogsr_All_2p))
Hog_All_3p_tt = unique(c(Hog_All_3p_t, Hogsr_All_3p))
Hog_All_4p_tt = unique(c(Hog_All_4p_t, Hogsr_All_4p))
Hog_All_5p_tt = unique(c(Hog_All_5p_t, Hogsr_All_5p))
Hog_All_6p_tt = unique(c(Hog_All_6p_t, Hogsr_All_6p))
Hog_All_7p_tt = unique(c(Hog_All_7p_t, Hogsr_All_7p))



do_full_norm_pca(subset(data_allso,subset = hog %in% Hog_All_1p_tt,
                   select = colnames(data_allso) %ni% c("Kf_Re2","Nc_Re2")),
            trait = "caste")
dim(subset(data_allso,subset = hog %in% Hog_All_1p_tt))
do_full_norm_pca(subset(data_allso,subset = hog %in% Hog_All_1p_tt,
                        select = colnames(data_allso) %ni% c("Kf_Re2","Nc_Re2")),
                 trait = "caste")


do_full_pca(subset(data_linso,subset = hog %in% Hogsr_Lin_1p,
                   select = colnames(data_linso) %ni% c("Kf_Re2","Nc_Re2")),
            trait = "caste")

do_full_norm_pca(subset(data_linso,subset = hog %in% Hogsr_Lin_1p,
                        select = colnames(data_linso) %ni% c()),
                 trait = "caste")

do_full_pca(subset(data_bifso,subset = hog %in% Hog_Bif_1p_tt,
                   select = colnames(data_bifso) %ni% c("Ap_So")),
            trait = "species")

do_full_norm_pca(subset(data_bifso,subset = hog %in% Hog_Bif_1p_tt,
                        select = colnames(data_bifso) %ni% c("Ap_So")
                        ),
                 trait = "caste")

##

saveRDS(c(data_allso,data_linso,data_bifso,
          Hog_All_1p_tt,Hog_All_2p_tt,Hog_All_3p_tt,
          Hog_All_4p_tt,Hog_All_5p_tt,Hog_All_6p_tt,
          Hog_All_7p_tt,Hogsr_Lin_4p,Hogsr_Lin_3p,
          Hogsr_Lin_2p,Hogsr_Lin_1p,Hog_Bif_5p_tt,
          Hog_Bif_4p_tt,Hog_Bif_3p_tt,Hog_Bif_2p_tt,
          Hog_Bif_1p_tt),file = "dataset_for_termite_PCA.rds")

pdf("norm_PCA_All_WoRe_9sp_1p_shared_genes_species.pdf")
do_full_norm_pca(subset(data_linsoold,subset = hog %in% Hogsr_Lin_1p,
                        select = colnames(data_linso) %ni% c()),
                 trait = "caste")
do_full_norm_pca(subset(data_bifso,subset = hog %in% Hog_Bif_1p_tt,
                        select = colnames(data_bifso) %ni% c("Ap_So")),
                 trait = "caste")
do_full_norm_pca(subset(data_allso,subset = hog %in% Hog_All_1p_tt,
                        select = colnames(data_allso) %ni% c("Ap_So")),
                 trait = "caste")
dim(subset(data_allso,subset = hog %in% Hog_All_1p_tt,
select = colnames(data_allso) %ni% c("Ap_So")))

dev.off()


pdf("norm_PCA_All_WoReSo_9sp_1p_shared_genes_caste.pdf")
do_full_norm_pca(subset(data_allso,subset = hog %in% Hog_All_1p_tt,
                        select = colnames(data_allso) %ni% c("Ap_So")),
                 trait = "caste")
dev.off()

pdf("norm_PCA_All_WoReSo_9sp_1p_shared_genes_species.pdf")
do_full_norm_pca(subset(data_allso,subset = hog %in% Hog_All_1p_tt,
                        select = colnames(data_allso) %ni% c("Ap_So")),
                 trait = "species")
dev.off()


pdf("norm_PCA_All_WoReSo_9sp_2p_shared_genes_caste.pdf")
do_full_norm_pca(subset(data_allso,subset = hog %in% Hog_All_2p_tt,
                        select = colnames(data_allso) %ni% c("Ap_So")),
                 trait = "caste")
dev.off()


pdf("norm_PCA_All_WoReSo_9sp_3p_shared_genes_caste.pdf")
do_full_norm_pca(subset(data_allso,subset = hog %in% Hog_All_3p_tt,
                        select = colnames(data_allso) %ni% c("Ap_So")),
                 trait = "caste")
dev.off()

dims = c(
dim(subset(data_allso,subset = hog %in% Hog_All_1p_tt,
           select = colnames(data_allso) %ni% c("Ap_So")))[1],
dim(subset(data_allso,subset = hog %in% Hog_All_2p_tt,
           select = colnames(data_allso) %ni% c("Ap_So")))[1],
dim(subset(data_allso,subset = hog %in% Hog_All_3p_tt,
           select = colnames(data_allso) %ni% c("Ap_So")))[1],
dim(subset(data_allso,subset = hog %in% Hog_All_4p_tt,
           select = colnames(data_allso) %ni% c("Ap_So")))[1],
dim(subset(data_allso,subset = hog %in% Hog_All_5p_tt,
           select = colnames(data_allso) %ni% c("Ap_So")))[1])

dims
# 5496 3512 1582  622  238

#

pdf("norm_PCA_Lin_WoReSo_4sp_1p_shared_genes_caste.pdf")
do_full_norm_pca(subset(data_linso,subset = hog %in% Hogsr_Lin_1p,
                        select = colnames(data_linso) %ni% c("Ap_So")),
                 trait = "caste")
dev.off()

pdf("norm_PCA_Bif_WoReSo_5sp_1p_shared_genes_caste.pdf")
do_full_norm_pca(subset(data_bifso,subset = hog %in% Hog_Bif_1p_tt,
                        select = colnames(data_bifso) %ni% c("Ap_So")),
                 trait = "caste")
dev.off()

pdf("norm_PCA_Lin_WoReSo_4sp_1p_shared_genes_species.pdf")
do_full_norm_pca(subset(data_linso,subset = hog %in% Hogsr_Lin_1p,
                        select = colnames(data_linso) %ni% c("Ap_So")),
                 trait = "species")
dev.off()

pdf("norm_PCA_Bif_WoReSo_5sp_1p_shared_genes_species.pdf")
do_full_norm_pca(subset(data_bifso,subset = hog %in% Hog_Bif_1p_tt,
                        select = colnames(data_bifso) %ni% c("Ap_So")),
                 trait = "species")
dev.off()
dims = c(
  dim(subset(data_linso,subset = hog %in% Hogsr_Lin_1p,
             select = colnames(data_linso) %ni% c("Ap_So")))[1],
  dim(subset(data_bifso,subset = hog %in% Hog_Bif_1p_tt,
             select = colnames(data_bifso) %ni% c("Ap_So")))[1])
dims
# 6211 5037



# ------------------------------------------------------------------------------
# ------------------ SEX -------------------------------------------------------
# ------------------------------------------------------------------------------


set_list = c("All","Lin", "Bif","Comp")
for (data_set in set_list){
  if (data_set == "Comp"){
    fin_df = SOG_Comp_13_species
    species_list = species_list_complete
  }
  if (data_set == "All"){
    fin_df = SOG_All_9_termites
    species_list = species_list_All
  }
  if (data_set == "Lin"){
    fin_df = SOG_Lin_termites
    species_list = species_list_Lin
  }
  if (data_set == "Bif"){
    fin_df = SOG_Bif_termites
    species_list = species_list_Bif
  }
  for (species in species_list){
    # species = "Anoplotermes_pacificus"
    # species = "Neotermes_castaneus"
    nreadcount = read.table(paste0("data/normalised_count_afterBC/",
                                   species,
                                   ".nreadcount.txt"))
    
    nreadcount$gene = rownames(nreadcount)
    nreadcount = left_join(nreadcount, gene2hog, by="gene")
    if (data_set == "Comp"){
      sogcount = subset(nreadcount, hog %in% SOG_Comp_13_species$hog)
    }
    if (data_set == "All"){
      sogcount = subset(nreadcount, hog %in% SOG_All_9_termites$hog)
    }
    if (data_set == "Lin"){
      sogcount = subset(nreadcount, hog %in% SOG_Lin_termites$hog)
    }
    if (data_set == "Bif"){
      sogcount = subset(nreadcount, hog %in% SOG_Bif_termites$hog)
    }

      list_col = colnames(sogcount)[grep("Fe", colnames(sogcount))]
      sogcount$Fe = rowMeans(sogcount[,colnames(sogcount) %in% list_col])
      list_col = colnames(sogcount)[grep("Ma", colnames(sogcount))]
      sogcount$Ma = rowMeans(sogcount[,colnames(sogcount) %in% list_col])
      sp = unlist(strsplit(colnames(nreadcount)[1],"_"))[1]
      inter_df = subset(sogcount, select = c(hog, Fe, Ma))
      colnames(inter_df) = c("hog", paste0(sp,"_Fe"),paste0(sp,"_Ma"))

  fin_df = left_join(fin_df, inter_df, by="hog")
}
  saveRDS(fin_df, paste0("data/data_for_",data_set,"_species_FeMa.RDS"))
}




data_all = readRDS("data/data_for_All_species_FeMa.RDS")

pdf("norm_PCA_all_FeMa_8sp_0p_shared_genes_sex.pdf")
do_full_norm_pca_sex(data_all,
            trait = "caste")
dev.off()

pdf("norm_PCA_all_FeMa_8sp_0p_shared_genes_species.pdf")
do_full_norm_pca_sex(data_all,
                     trait = "species")
dev.off()

# here need to modify manually the function to remove Ap
pdf("norm_PCA_all_FeMa_7sp_0p_shared_genes_sex.pdf")
do_full_norm_pca_sex(data_all,
                     trait = "caste")
dev.off()

pdf("norm_PCA_all_FeMa_7sp_0p_shared_genes_species.pdf")
do_full_norm_pca_sex(data_all,
                     trait = "species")
dev.off()


do_full_norm_pca_sex = function(data,trait, dimA=1, dimB=2){
  data_scaled = data[1]
  spelist = unique(unlist(strsplit(colnames(data[-1]),"_",fixed = T)))
  spelist = spelist[spelist %ni% c("Wo","Re","Re2","So", "Fe", "Ma","Mn", "Ap")]
  for (sp in spelist){
    data_mean = subset(data,
                       select = c("hog",colnames(data)[grep(sp, colnames(data))]))
    row.names(data_mean) = data_mean$hog
    #M = mean(rowMeans(data_mean[-1]))
    #Ncol = dim(data_mean[-1])[2]
    #center_vector = rep(M,Ncol)
    M = (rowMeans(data_mean[-1]))
    center_vector = (M)
    data_scaled_temp = data.frame(scale(t(data_mean[-1]), center = center_vector, scale = center_vector))
    data_scaled_temp = data.frame(t(data_scaled_temp))
    data_scaled_temp$hog  =row.names(data_scaled_temp)
    data_scaled = left_join(data_scaled,data_scaled_temp, by= "hog")
    
  }
  
  
  
  data1=data_scaled[-1]
  tdata = t(data1)
  colnames(tdata) = data$hog
  species = t(data.frame(strsplit(row.names(tdata),split = "_"))[1,])
  caste = t(data.frame(strsplit(row.names(tdata),split = "_"))[2,])
  sex = t(data.frame(strsplit(row.names(tdata),split = "_"))[3,])
  colnames(species) = "species"
  colnames(caste) = "caste"
  colnames(sex) = "sex"
  row.names(species) = row.names(tdata)
  row.names(caste) = row.names(tdata)
  row.names(sex) = row.names(tdata)
  data = cbind.data.frame(tdata,species,caste,sex)
  title_pca=""
  res.pca = PCA(data[,-c(ncol(data)-2,ncol(data)-1,ncol(data))],
                graph = FALSE,
                scale.unit = F)
  if(trait == "species"){
    RES = fviz_pca_ind(res.pca,  axes = c(dimA,dimB),
                       label="none", habillage=as.factor(data$species),
                       addEllipses=TRUE, ellipse.level=0.7, title=title_pca)
  }
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
  RES
}


###################
#END#
###################


