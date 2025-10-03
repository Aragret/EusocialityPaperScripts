################################################################################
#
#             Title: SET for Caste Specific Genes 
#           Project: Eusociality
#           Authors: C. Aumont
#              Year: 2025
#
################################################################################


#-------------------------------------------------------------------------------
#----------- Introduction ------------------------------------------------------
#-------------------------------------------------------------------------------

#----------- Library installation -----------------------------------------------

library("RColorBrewer")
library(readxl)
library(SuperExactTest)
library(dplyr)
"%ni%" = Negate("%in%")

# ---------- Set work directory ------------------------------------------------

# enter the working directory. It must include the data folder.
setwd("~/WORK/Dino_McMahon/Work/Lab work and projects/6. Gene network/10.SET/2.CSG")


# ---------- Uploading annex files -------------------------------------------

#------- gene2hog --------------------------------------------------------------

gene2hog = read_xlsx(path = "../../TABLES/Tables_SA.xlsx",
                     sheet = "S2_gene_to_hog",
                     skip = 2)

#------- OR ID -----------------------------------------------------------------

OG2gene = read.table("../../ORs/OR/gene2OG_map.txt",
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

# formating HOG column
OG2gene$hog = gsub("OG", "OrOG", OG2gene$orthogroup)

summary(OG2gene)
head(OG2gene)


#------- gene to hog afterBC ---------------------------------------------------


# adding Genes that are not in gene2hog
subOG2gene = subset(OG2gene, spe %in% c("Apac","Cges",
                                        "Hsjo","Kfla", 
                                        "Mdar","Mnat",
                                        "Ncas","PRsim",
                                        "Rfla","Znev",
                                        "Cmer","Cpun","BGER","Bori"))

genetoadd = subOG2gene[subOG2gene$Node %ni% gene2hog$gene,c(4,5)]
colnames(genetoadd) = c("gene", "hog")
gene2hog = bind_rows(gene2hog, genetoadd)

# replacing HOG ids by OrHOG ids in gene to hog file
for (gene in gene2hog$gene){
  if (gene %in% OG2gene$Node){
    gene2hog$hog[gene2hog$gene==gene] = OG2gene$hog[OG2gene$Node == gene] 
  }
}
tail(gene2hog)
write.table(gene2hog, "gene2hog_afterBC.txt",
            sep="\t",
            quote = F,
            row.names = F)

# TODO: add cluster2 genes to the list

OGclus2 = read.table("../../ORs/OR/Cluster2_genenames.txt")
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

# replacing HOG ids by OrHOG ids in gene to hog file
for (gene in gene2hog$gene){
  if (gene %in% OGclus2$gene){
    gene2hog$hog[gene2hog$gene==gene] = "cluster2" 
  }
}
tail(gene2hog,50)
write.table(gene2hog, "gene2hog_afterBC_inclClus2.txt",
            sep="\t",
            quote = F,
            row.names = F)




#------- Master sheet ----------------------------------------------------------

master = data.frame(read_excel("../../TABLES/small_master_sheet_eusoc.xlsx",
                               sheet = "S0_Master_sheet",
                               na = c("","NA")))

# list of single copy orthologues for the 9 termite species
SOG_list_9_termites = data.frame(
  hog=read.table("../1.DEG/Pairwise_WORE/Script/SOG_list_9spe.txt")$V1,
  sog="yes")



spe_list = c("Apac","Cges","Hsjo","Kfla", "Mdar","Mnat","Ncas","PRsim","Rfla",
             "Znev","Cmer","Cpun","BGER","Bori")
abg_list = c("g","a")
caste_list = c("Wo","So","Ny","Al","Pr","Er","Nm","Ne","Ju","Ad")
species_list = c(               # Attention: the order matters: need to be the same as spe_list
  "Anoplotermes_pacificus",
  "Coptotermes_gestroi",
  "Hodotermopsis_sjostedti",
  "Kalotermes_flavicollis",
  "Mastotermes_darwiniensis",
  "Macrotermes_natalensis",
  "Neotermes_castaneus",
  "Prorhinotermes_simplex",
  "Reticulitermes_flavipes",
  "Zootermopsis_nevadensis",
  "Cryptocercus_meridianus",
  "Cryptocercus_punctulatus",
  "Blattella_germanica",
  "Blatta_orientalis")



# ------------------------------------------------------------------------------
# ---------- functions ---------------------------------------------------------
# ------------------------------------------------------------------------------

# Return the list of shared hog among the number of species asked (degree)
return_common_hogs = function(res, degree){
  #  Result=supertest(data, n)
  #  res = summary(Result)$Table
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
  return(unique(hog_list))
}

# Return the list of shared hogs with Mastotermes included
return_common_hogs_3M = function(res, degree){
  #  Result=supertest(data, n)
  #  res = summary(Result)$Table
  deg_res =res[res$Degree >= degree,]
  deg_resM = deg_res$Elements[deg_res$Intersections %in% deg_res$Intersections[grepl("Mdar",deg_res$Intersections)]]
  hog_list =c()
  for (elem in deg_resM){
    hog_line = unlist(strsplit(elem,", "))
    for (hog in hog_line){
      if (hog %ni% hog_list){
        hog_list = append(hog_list, hog)
      }  
    }
  }
  return(unique(hog_list))
}

# ------------------------------------------------------------------------------
# ---------- n calculation ---------------------------------------------------------
# ------------------------------------------------------------------------------

# here n is the number of shared hogs between the species included in the comparison

# Solitary cockroaches
hog_BGER = unique(subset(gene2hog,subset = grepl("BGER",x = gene2hog$gene))$hog)
hog_Bori = unique(subset(gene2hog,subset = grepl("Bori",x = gene2hog$gene))$hog)
n_soli_ck = length(SuperExactTest::intersect(hog_BGER,hog_Bori))
soli_ck = SuperExactTest::intersect(hog_BGER,hog_Bori)

# subsocial cockroaches
hog_Cmer = unique(subset(gene2hog,subset = grepl("Cmer",x = gene2hog$gene))$hog)
hog_Cpun = unique(subset(gene2hog,subset = grepl("Cpun",x = gene2hog$gene))$hog)
n_sub_ck = length(SuperExactTest::intersect(hog_Cmer,hog_Cmer))
sub_ck = SuperExactTest::intersect(hog_Cmer,hog_Cmer)

# Termites
hog_Mnat = unique(subset(gene2hog,subset = grepl("Mnat",x = gene2hog$gene))$hog)
hog_Mdar = unique(subset(gene2hog,subset = grepl("Mdar",x = gene2hog$gene))$hog)
hog_Apac = unique(subset(gene2hog,subset = grepl("Apac",x = gene2hog$gene))$hog)
hog_Ncas = unique(subset(gene2hog,subset = grepl("Ncas",x = gene2hog$gene))$hog)
hog_Kfla = unique(subset(gene2hog,subset = grepl("Kfla",x = gene2hog$gene))$hog)
hog_Rfla = unique(subset(gene2hog,subset = grepl("Rfla",x = gene2hog$gene))$hog)
hog_Cges = unique(subset(gene2hog,subset = grepl("Cges",x = gene2hog$gene))$hog)
hog_Hsjo = unique(subset(gene2hog,subset = grepl("Hsjo",x = gene2hog$gene))$hog)
hog_PRsim = unique(subset(gene2hog,subset = grepl("PRsim",x = gene2hog$gene))$hog)
hog_Znev = unique(subset(gene2hog,subset = grepl("Znev",x = gene2hog$gene))$hog)

linter=  SuperExactTest::intersect(hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim, hog_Znev)
n_linter = length(linter)

bifter = SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Apac,hog_Cges,hog_Rfla)
n_bifter = length(bifter)
                 

x10ter = SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Apac,hog_Cges,hog_Rfla,
                                          hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim,hog_Znev)
n10ter = length(x10ter)
x9ter = SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Cges,hog_Rfla,
                                          hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim,hog_Znev)
n9ter = length(x9ter)
Nyter = SuperExactTest::intersect(hog_Rfla,hog_Kfla,hog_Ncas,hog_Znev)
nNyter = length(Nyter)
Prter = SuperExactTest::intersect(hog_Cges,hog_PRsim,hog_Apac)
nPrter = length(Prter)
Srter = SuperExactTest::intersect(hog_Mdar,hog_Rfla,hog_Hsjo,hog_Kfla,hog_Ncas)
nSrter = length(Srter)
Alter = SuperExactTest::intersect(hog_Mnat,hog_Kfla,hog_Ncas)
nAlter = length(Alter)


# ------------------------------------------------------------------------------
# ---------- main ---------------------------------------------------------
# ------------------------------------------------------------------------------
datas = list()

for (caste in caste_list){

  for (spe in spe_list){
    species = species_list[match(spe,spe_list)]
    file_name = paste0("data/Node_groups_normalised_afterBC/Node_groups_normalisedAfterBC_",
                      species,".csv")
    file_of_interest = read.csv2(file = file_name)
    abg = "g" 
    group_of_interest = paste0(abg,".Network_",species,"_", caste)
    gene_of_interest = file_of_interest$Node[
      file_of_interest$Group_pval_Phi_tilde == group_of_interest
      ]
    if (length(gene_of_interest) > 0){
     # df = data.frame(spec= rep(spe,length(gene_of_interest)), goi=gene_of_interest)
    #  df$HOG = NA
    #  for (gene in df$goi){
    #    if (gene %in% gene2hog$gene){
    #      df$HOG[df$goi == gene] = gene2hog$hog[gene2hog$gene == gene]
    #    }else{
    #      df$HOG[df$goi == gene] = NA
    #    }
    #  }
      ###
      df = data.frame(spec= rep(spe,length(gene_of_interest)), gene=gene_of_interest)
      df1 = left_join(df, gene2hog, by= "gene")
      Temp =list(V1=unique(na.omit(df1$hog)))
      ###
      #Temp =list(V1=unique(na.omit(df$HOG)))
      names(Temp) = paste0(spe, "_", caste)
      datas = append(datas, Temp)
    }
  }
  print(names(datas))   

} 

saveRDS(datas, file= "Result/datas_SET_CSG_afterBC_withclus2.rds")
datas = readRDS("Result/datas_SET_CSG_afterBC_withclust2.rds")

# ------------------------------------------------------------------------------
# ------- SET and shared hog extraction ----------------------------------------
# ------------------------------------------------------------------------------


# ---------- Final table preparation -------------------------------------------
Final_SET_Table = data.frame(caste=character(),
                             Intersections=character(),
                             Degree = integer(),
                             Observed.Overlap = numeric(),
                             Expected.Overlap = numeric(),
                             FE = numeric(),
                             P.value = numeric(),
                             Elements = character())

Final_SharedHOG_table = data.frame(hog=character(),
                                   type=character()
)



statistic = F
# ------- All worker termites ----------------------------------------

data_Wo_termites = datas[names(datas) %in% c(
  "Apac_Wo","Cges_Wo","Hsjo_Wo","Kfla_Wo", "Mdar_Wo","Mnat_Wo",
  "Ncas_Wo","PRsim_Wo","Rfla_Wo","Znev_Wo")]
if (statistic) {
  dstat = data_Wo_termites
  data_Wo_termites = list(
    Apac_Wo = subset(dstat$Apac_Wo,dstat$Apac_Wo %in% x10ter),
    Cges_Wo = subset(dstat$Cges_Wo,dstat$Cges_Wo %in% x10ter),
    Hsjo_Wo = subset(dstat$Hsjo_Wo,dstat$Hsjo_Wo %in% x10ter),
    Kfla_Wo = subset(dstat$Kfla_Wo,dstat$Kfla_Wo %in% x10ter),
    Mdar_Wo = subset(dstat$Mdar_Wo,dstat$Mdar_Wo %in% x10ter),
    Mnat_Wo = subset(dstat$Mnat_Wo,dstat$Mnat_Wo %in% x10ter),
    Ncas_Wo = subset(dstat$Ncas_Wo,dstat$Ncas_Wo %in% x10ter),
    PRsim_Wo = subset(dstat$PRsim_Wo,dstat$PRsim_Wo %in% x10ter),
    Znev_Wo = subset(dstat$Znev_Wo,dstat$Znev_Wo %in% x10ter),
    Rfla_Wo = subset(dstat$Rfla_Wo,dstat$Rfla_Wo %in% x10ter)
  )
  Result=supertest(data_Wo_termites, n = n10ter)
} else {
  Result=supertest(data_Wo_termites)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Wo_allter"
Final_SET_Table = bind_rows(Final_SET_Table, res)
#rch = data.frame(type = "Wo_10p",
 #                hog = return_common_hogs(res = res,degree = 10))
#Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Wo_9p",
                 hog = return_common_hogs(res = res,degree = 9))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Wo_7p",
                 hog = return_common_hogs(res = res,degree = 7))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Wo_5p",
                 hog = return_common_hogs(res = res,degree = 5))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


data_So_termites = datas[names(datas) %in% c(
  "Cges_So","Hsjo_So","Kfla_So", "Mdar_So","Mnat_So",
  "Ncas_So","PRsim_So","Rfla_So","Znev_So")]
if (statistic) {
  dstat = data_So_termites
  data_So_termites = list(
    Cges_So = subset(dstat$Cges_So,dstat$Cges_So %in% x9ter),
    Hsjo_So = subset(dstat$Hsjo_So,dstat$Hsjo_So %in% x9ter),
    Kfla_So = subset(dstat$Kfla_So,dstat$Kfla_So %in% x9ter),
    Mdar_So = subset(dstat$Mdar_So,dstat$Mdar_So %in% x9ter),
    Mnat_So = subset(dstat$Mnat_So,dstat$Mnat_So %in% x9ter),
    Ncas_So = subset(dstat$Ncas_So,dstat$Ncas_So %in% x9ter),
    PRsim_So = subset(dstat$PRsim_So,dstat$PRsim_So %in% x9ter),
    Znev_So = subset(dstat$Znev_So,dstat$Znev_So %in% x9ter),
    Rfla_So = subset(dstat$Rfla_So,dstat$Rfla_So %in% x9ter)
  )
  Result=supertest(data_So_termites, n = n9ter)
} else {
  Result=supertest(data_So_termites)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="So_allter"
Final_SET_Table = bind_rows(Final_SET_Table, res)
#rch = data.frame(type = "So_9p",
 #                hog = return_common_hogs(res = res,degree = 9))
#Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "So_7p",
                 hog = return_common_hogs(res = res,degree = 7))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "So_5p",
                 hog = return_common_hogs(res = res,degree = 5))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)





data_Ny_termites = datas[names(datas) %in% c(
 "Kfla_Ny", "Ncas_Ny","Rfla_Ny","Znev_Ny")]
if (statistic) {
  dstat = data_Ny_termites
  data_Ny_termites = list(
    Kfla_Ny = subset(dstat$Kfla_Ny,dstat$Kfla_Ny %in% Nyter),
    Ncas_Ny = subset(dstat$Ncas_Ny,dstat$Ncas_Ny %in% Nyter),
    Znev_Ny = subset(dstat$Znev_Ny,dstat$Znev_Ny %in% Nyter),
    Rfla_Ny = subset(dstat$Rfla_Ny,dstat$Rfla_Ny %in% Nyter)
  )
  Result=supertest(data_Ny_termites, n = nNyter)
} else {
  Result=supertest(data_Ny_termites)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Ny_allter"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "Ny_4p",
                 hog = return_common_hogs(res = res,degree = 4))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Ny_3p",
                 hog = return_common_hogs(res = res,degree = 3))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)





data_Pr_termites = datas[names(datas) %in% c(
  "Apac_Pr", "Cges_Pr","PRsim_Pr")]
if (statistic) {
  dstat = data_Pr_termites
  data_Pr_termites = list(
    Apac_Pr = subset(dstat$Apac_Pr,dstat$Apac_Pr %in% Prter),
    Cges_Pr = subset(dstat$Cges_Pr,dstat$Cges_Pr %in% Prter),
    PRsim_Pr = subset(dstat$PRsim_Pr,dstat$PRsim_Pr %in% Prter)
  )
  Result=supertest(data_Pr_termites, n = nPrter)
} else {
  Result=supertest(data_Pr_termites)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Pr_allter"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "Pr_3p",
                 hog = return_common_hogs(res = res,degree = 3))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)





data_Al_termites = datas[names(datas) %in% c(
  "Kfla_Al","Mnat_Al", "Ncas_Al")]
if (statistic) {
  dstat = data_Al_termites
  data_Al_termites = list(
    Kfla_Al = subset(dstat$Kfla_Al,dstat$Kfla_Al %in% Alter),
    Mnat_Al = subset(dstat$Mnat_Al,dstat$Mnat_Al %in% Alter),
    Ncas_Al = subset(dstat$Ncas_Al,dstat$Ncas_Al %in% Alter)
  )
  Result=supertest(data_Al_termites, n = nAlter)
} else {
  Result=supertest(data_Al_termites)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Al_allter"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "Al_3p",
                 hog = return_common_hogs(res = res,degree = 3))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


data_Sr_termites = datas[names(datas) %in% c(
  "Hsjo_Er","Kfla_Er", "Mdar_Er","Ncas_Er","Rfla_Nm")]
if (statistic) {
  dstat = data_Sr_termites
  data_Sr_termites = list(
    Hsjo_Er = subset(dstat$Hsjo_Er,dstat$Hsjo_Er %in% Srter),
    Kfla_Er = subset(dstat$Kfla_Er,dstat$Kfla_Er %in% Srter),
    Mdar_Er = subset(dstat$Mdar_Er,dstat$Mdar_Er %in% Srter),
    Ncas_Er = subset(dstat$Ncas_Er,dstat$Ncas_Er %in% Srter),
    Rfla_Nm = subset(dstat$Rfla_Nm,dstat$Rfla_Nm %in% Srter)
  )
  Result=supertest(data_Sr_termites, n = nSrter)
} else {
  Result=supertest(data_Sr_termites)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Sr_allter"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "Sr_4p",
                 hog = return_common_hogs(res = res,degree = 4))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Sr_3p",
                 hog = return_common_hogs(res = res,degree = 3))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)




data_Ju_Subcockroaches = datas[names(datas) %in% c(
  "Cmer_Ju","Cpun_Ju")]
if (statistic) {
  dstat = data_Ju_Subcockroaches
  data_Ju_Subcockroaches = list(
    Cmer_Ju = subset(dstat$Cmer_Ju,dstat$Cmer_Ju %in% sub_ck),
    Cpun_Ju = subset(dstat$Cpun_Ju,dstat$Cpun_Ju %in% sub_ck)
  )
  Result=supertest(data_Ju_Subcockroaches, n = n_sub_ck)
} else {
  Result=supertest(data_Ju_Subcockroaches)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Ju_Subck"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "JuSu_2p",
                 hog = return_common_hogs(res = res,degree = 2))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)

data_Ju_Solcockroaches = datas[names(datas) %in% c(
  "BGER_Ju","Bori_Ju")]
if (statistic) {
  dstat = data_Ju_Solcockroaches
  data_Ju_Solcockroaches = list(
    BGER_Ju = subset(dstat$BGER_Ju,dstat$BGER_Ju %in% soli_ck),
    Bori_Ju = subset(dstat$Bori_Ju,dstat$Bori_Ju %in% soli_ck)
  )
  Result=supertest(data_Ju_Solcockroaches, n = n_soli_ck)
} else {
  Result=supertest(data_Ju_Solcockroaches)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Ju_Solck"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "JuSo_2p",
                 hog = return_common_hogs(res = res,degree = 2))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


data_Ad_Subcockroaches = datas[names(datas) %in% c(
  "Cmer_Ad","Cpun_Ad")]
if (statistic) {
  dstat = data_Ad_Subcockroaches
  data_Ad_Subcockroaches = list(
    Cmer_Ad = subset(dstat$Cmer_Ad,dstat$Cmer_Ad %in% sub_ck),
    Cpun_Ad = subset(dstat$Cpun_Ad,dstat$Cpun_Ad %in% sub_ck)
  )
  Result=supertest(data_Ad_Subcockroaches, n = n_sub_ck)
} else {
  Result=supertest(data_Ad_Subcockroaches)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Ad_Subck"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "AdSu_2p",
                 hog = return_common_hogs(res = res,degree = 2))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)

data_Ad_Solcockroaches = datas[names(datas) %in% c(
  "BGER_Ad","Bori_Ad")]
if (statistic) {
  dstat = data_Ad_Solcockroaches
  data_Ad_Solcockroaches = list(
    BGER_Ad = subset(dstat$BGER_Ad,dstat$BGER_Ad %in% soli_ck),
    Bori_Ad = subset(dstat$Bori_Ad,dstat$Bori_Ad %in% soli_ck)
  )
  Result=supertest(data_Ad_Solcockroaches, n = n_soli_ck)
} else {
  Result=supertest(data_Ad_Solcockroaches)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Ad_Solck"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "AdSo_2p",
                 hog = return_common_hogs(res = res,degree = 2))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------------------------------------------------------------------------------
# ----------- Create shared hog table ------------------------------------------
# ------------------------------------------------------------------------------

# ----------- extract table information ----------------------------------------

master$hog = as.factor(master$hog)

# ----------- Merge tables -----------------------------------------------------
spelist = c()
V1 = data.frame(do.call('rbind', strsplit(as.character(gene2hog$gene), '0', fixed = T)))[,1]
V2 = data.frame(do.call('rbind', strsplit(as.character(V1), 'O', fixed = T)))[,1]
spelist = gsub("Bger", "BGER", V2)
gene2hog$spe = spelist

B = left_join(x = Final_SharedHOG_table,
              y = subset(gene2hog, spe %ni% c("Nluj")),
              by="hog")
master$hog <-NULL
C = left_join(x = B, y = master, by="gene")


D = left_join(x = C, y = SOG_list_9_termites, by="hog")


# ------------------------------------------------------------------------------
# ----------- Saving files -----------------------------------------------------
# ------------------------------------------------------------------------------



# ------------ Save the super exact test results -------------------------------

if (statistic){
  write.table(Final_SET_Table, 
              "Result/SuperExactTest_caste-specific-genes_AfterBC_stats_withclus2.txt", 
              quote = F,
              row.names = F,
              sep = "\t",
              dec = ",")
}else{
  write.table(Final_SET_Table,
              "Result/SuperExactTest_caste-specific-genes_AfterBC_nostats_withclus2.txt", 
              quote = F,
              row.names = F,
              sep = "\t",
              dec = ",")
}

# ------------ Save the shared hog results -------------------------------------
if (!statistic){
  write.table(D,file = "Result/Shared_HOG_caste-specific-genes_AfterBC_nostat_withclus2.txt",
              quote = F,sep="\t",
              row.names = F,col.names = T)
}
if (statistic){
  write.table(D,file = "Result/Shared_HOG_caste-specific-genes_AfterBC_stat_withclus2.txt",
              quote = F,sep="\t",
              row.names = F,col.names = T)
}

