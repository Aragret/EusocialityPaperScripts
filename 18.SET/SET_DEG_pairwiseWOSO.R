################################################################################
#
#             Title: SET for DEG in Worker-Reproductive pairwise comparison 
#           Project: Eusociality
#           Authors: C. Aumont
#              Year: 2025-07
#
################################################################################

#-------------------------------------------------------------------------------
#----------- Introduction ------------------------------------------------------
#-------------------------------------------------------------------------------

#----------- Library installation ----------------------------------------------

library("RColorBrewer")
library(readxl)
library(dplyr)
library(SuperExactTest)
"%ni%" = Negate("%in%")


# ---------- Set work directory ------------------------------------------------

# enter the working directory. it must include the folder with the Node_groups
setwd("~/WORK/Dino_McMahon/Work/Lab work and projects/6. Gene network/10.SET/1.DEG/Pairwise_WOSO/Script/")


# ---------- Final table preparation -------------------------------------------

datas_Wo = list()
datas_So = list()

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


# ---------- Uploading annex files -------------------------------------------

# gene ID to HOG ID file
gene2hog = read.table("../../../2.CSG/gene2hog_afterBC.txt",
                      header = T)


# DEG pairwise category table
data_raw = data.frame(read_xlsx("../../../../TABLES/Tables_SA.xlsx",
                                sheet = "DEG_pairwise&specif_aftBC_0.05",
                                col_types = "text",na = "NA"))
data_raw$Species = factor(data_raw$Species)
data_raw$Al_Ne = factor(data_raw$Al_Ne)
data_raw$Al_Ny = factor(data_raw$Al_Ny)
data_raw$Al_So = factor(data_raw$Al_So)
data_raw$Al_Wo = factor(data_raw$Al_Wo)
data_raw$Ne_Ny = factor(data_raw$Ne_Ny)
data_raw$Ne_So = factor(data_raw$Ne_So)
data_raw$Ne_Wo = factor(data_raw$Ne_Wo)
data_raw$Ny_So = factor(data_raw$Ny_So)
data_raw$Ny_Wo = factor(data_raw$Ny_Wo)
data_raw$Er_So = factor(data_raw$Er_So)
data_raw$So_Wo = factor(data_raw$So_Wo)
data_raw$Er_Wo = factor(data_raw$Er_Wo)
data_raw$Nm_Ny = factor(data_raw$Nm_Ny)
data_raw$Nm_So = factor(data_raw$Nm_So)
data_raw$Nm_Wo = factor(data_raw$Nm_Wo)
data_raw$Pr_So = factor(data_raw$Pr_So)
data_raw$Pr_Wo = factor(data_raw$Pr_Wo)
data_raw$Ad_Ju = factor(data_raw$Ad_Ju)

# Master sheet
master = data.frame(read_excel("../../../../TABLES/small_master_sheet_eusoc.xlsx",
                               sheet = "S0_Master_sheet",
                               na = c("","NA")))


# list of single copy orthologues for the 9 termite species
SOG_list_9_termites = data.frame(
  hog=read.table("SOG_list_9spe.txt")$V1,
  sog="yes")


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

# Termites
hog_Mnat = unique(subset(gene2hog,subset = grepl("Mnat",x = gene2hog$gene))$hog)
hog_Mdar = unique(subset(gene2hog,subset = grepl("Mdar",x = gene2hog$gene))$hog)
hog_Ncas = unique(subset(gene2hog,subset = grepl("Ncas",x = gene2hog$gene))$hog)
hog_Kfla = unique(subset(gene2hog,subset = grepl("Kfla",x = gene2hog$gene))$hog)
hog_Rfla = unique(subset(gene2hog,subset = grepl("Rfla",x = gene2hog$gene))$hog)
hog_Cges = unique(subset(gene2hog,subset = grepl("Cges",x = gene2hog$gene))$hog)
hog_Hsjo = unique(subset(gene2hog,subset = grepl("Hsjo",x = gene2hog$gene))$hog)
hog_PRsim = unique(subset(gene2hog,subset = grepl("PRsim",x = gene2hog$gene))$hog)
hog_Znev = unique(subset(gene2hog,subset = grepl("Znev",x = gene2hog$gene))$hog)

nlinter = length(SuperExactTest::intersect(hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim,hog_Znev))
nbifter = length(SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Cges,hog_Rfla))
nallter = length(SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Znev,hog_Cges,hog_Rfla,
                            hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim))
linter = SuperExactTest::intersect(hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim,hog_Znev)
bifter = SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Cges,hog_Rfla)
allter = SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Znev,hog_Cges,hog_Rfla,
                                            hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim)

# ---------------Initialisation ------------------------------------------------

spe_list = c("Cges","Hsjo","Kfla", "Mdar","Mnat","Ncas","PRsim","Rfla",
             "Znev")

species_list = c(               # Attention: the order matters: need to be the same as spe_list
  "Coptotermes_gestroi",
  "Hodotermopsis_sjostedti",
  "Kalotermes_flavicollis",
  "Mastotermes_darwiniensis",
  "Macrotermes_natalensis",
  "Neotermes_castaneus",
  "Prorhinotermes_simplex",
  "Reticulitermes_flavipes",
  "Zootermopsis_nevadensis")

colu_woso = c("So_Wo")

# ------------------------------------------------------------------------------

for (spe in spe_list){
  species = species_list[match(spe,spe_list)]
  print(spe)
  data_sub = subset(data_raw, subset= data_raw$Species == spe)
  sub_sub_data = subset(data_sub, select = colnames(data_sub) %in% c("gene_name",colu_woso))
  for (colu in colu_woso){
    if (!is.na(sub_sub_data[1,colu])){
      castes = unlist(strsplit(colu,split = "_"))
      for (casteA in castes){      
        df0 = subset(sub_sub_data, select = colnames(sub_sub_data) %in% c("gene_name",colu))
        list_caste_A = df0$gene_name[df0[2]==casteA]
        if (!is.na(list_caste_A[1])) {
          gene_of_interest = list_caste_A
        }
        if (length(gene_of_interest) > 0){
          df = data.frame(spec= rep(spe,length(gene_of_interest)),
                    gene=gene_of_interest)
        }
        df1 = left_join(df, gene2hog, by= "gene")

        Temp =list(V1=unique(na.omit(df1$hog)))
        if (casteA %in% c("Wo")){
          names(Temp) = paste0(spe,"_Wo")
          datas_Wo = append(datas_Wo, Temp)
        }
        if (casteA %ni% c("Wo")){
          names(Temp) = paste0(spe,"_So")
          datas_So = append(datas_So, Temp)
        }
      }
    }       
  }
}

saveRDS(datas_So, "../Result/data_biased_So.rds")
saveRDS(datas_Wo, "../Result/data_biased_Wo.rds")
names(datas_So)
# ------------------------------------------------------------------------------
# ------- SET and shared hog extraction ----------------------------------------
# ------------------------------------------------------------------------------
datas_Wo = readRDS("../Result/data_biased_Wo.rds")
datas_So = readRDS("../Result/data_biased_So.rds")
statistic = F

data_Wo_termites = datas_Wo
if (statistic) {
  dstat = data_Wo_termites
  data_Wo_termites = list(
    Cges_Wo = subset(dstat$Cges_Wo,dstat$Cges_Wo %in% allter),
    Hsjo_Wo = subset(dstat$Hsjo_Wo,dstat$Hsjo_Wo %in% allter),
    Kfla_Wo = subset(dstat$Kfla_Wo,dstat$Kfla_Wo %in% allter),
    Mdar_Wo = subset(dstat$Mdar_Wo,dstat$Mdar_Wo %in% allter),
    Mnat_Wo = subset(dstat$Mnat_Wo,dstat$Mnat_Wo %in% allter),
    Ncas_Wo = subset(dstat$Ncas_Wo,dstat$Ncas_Wo %in% allter),
    PRsim_Wo = subset(dstat$PRsim_Wo,dstat$PRsim_Wo %in% allter),
    Rfla_Wo = subset(dstat$Rfla_Wo,dstat$Rfla_Wo %in% allter),
    Znev_Wo = subset(dstat$Znev_Wo,dstat$Znev_Wo %in% allter)
  )
  Result=supertest(data_Wo_termites, n = nallter)
} else {
  Result=supertest(data_Wo_termites)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Wo_allter"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "Wo_9p",
                 hog = return_common_hogs(res = res,degree = 9))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Wo_7p",
                 hog = return_common_hogs(res = res,degree = 7))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Wo_5p",
                 hog = return_common_hogs(res = res,degree = 5))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- All Soldier termites ----------------------------------------

data_So_termites = datas_So
if (statistic) {
  dstat = data_So_termites
  data_So_termites = list(
    Cges_So = subset(dstat$Cges_So,dstat$Cges_So %in% allter),
    Hsjo_So = subset(dstat$Hsjo_So,dstat$Hsjo_So %in% allter),
    Kfla_So = subset(dstat$Kfla_So,dstat$Kfla_So %in% allter),
    Mdar_So = subset(dstat$Mdar_So,dstat$Mdar_So %in% allter),
    Mnat_So = subset(dstat$Mnat_So,dstat$Mnat_So %in% allter),
    Ncas_So = subset(dstat$Ncas_So,dstat$Ncas_So %in% allter),
    PRsim_So = subset(dstat$PRsim_So,dstat$PRsim_So %in% allter),
    Rfla_So = subset(dstat$Rfla_So,dstat$Rfla_So %in% allter),
    Znev_So = subset(dstat$Znev_So,dstat$Znev_So %in% allter)
  )
  Result=supertest(data_So_termites, n = nallter)
} else {
  Result=supertest(data_So_termites)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="So_allter"
Final_SET_Table = bind_rows(Final_SET_Table, res)
#rch = data.frame(type = "So_9p",
 #             hog = return_common_hogs(res = res,degree = 9)) 
#Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
#rch = data.frame(type = "So_7p",
 #                hog = return_common_hogs(res = res,degree = 7))
#Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "So_5p",
                 hog = return_common_hogs(res = res,degree = 5))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Worker termites with linear ontogeny ---------------------------

data_Wo_linear = datas_Wo[names(datas_Wo) %in% c("Hsjo_Wo", "Kfla_Wo",
                                                 "Ncas_Wo", "PRsim_Wo","Znev_Wo")]
if (statistic) {
  dstat = data_Wo_linear
  data_Wo_linear = list(
    Hsjo_Wo = subset(dstat$Hsjo_Wo,dstat$Hsjo_Wo %in% linter),
    Kfla_Wo = subset(dstat$Kfla_Wo,dstat$Kfla_Wo %in% linter),
    Ncas_Wo = subset(dstat$Ncas_Wo,dstat$Ncas_Wo %in% linter),
    PRsim_Wo = subset(dstat$PRsim_Wo,dstat$PRsim_Wo %in% linter),
    Znev_Wo = subset(dstat$Znev_Wo,dstat$Znev_Wo %in% linter)
  )
  Result=supertest(data_Wo_linear, n = nlinter)
} else {
  Result=supertest(data_Wo_linear)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Wo_linear"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "WoLin_5p",
                 hog = return_common_hogs(res = res,degree = 5))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "WoLin_4p",
                 hog = return_common_hogs(res = res,degree = 4))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)

# ------- Worker termites with bifurcated ontogeny -----------------------

data_Wo_bifur = datas_Wo[names(datas_Wo) %in% c("Mdar_Wo", "Rfla_Wo",
                                                "Cges_Wo", "Mnat_Wo")]
if (statistic) {
  dstat = data_Wo_bifur
  data_Wo_bifur = list(
    Cges_Wo = subset(dstat$Cges_Wo,dstat$Cges_Wo %in% bifter),
    Mdar_Wo = subset(dstat$Mdar_Wo,dstat$Mdar_Wo %in% bifter),
    Mnat_Wo = subset(dstat$Mnat_Wo,dstat$Mnat_Wo %in% bifter),
    Rfla_Wo = subset(dstat$Rfla_Wo,dstat$Rfla_Wo %in% bifter)
  )
  Result=supertest(data_Wo_bifur, n = nbifter)
} else {
  Result=supertest(data_Wo_bifur)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Wo_bifur"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "WoBif_4p",
                 hog = return_common_hogs(res = res,degree = 4))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)



# ------- Soldier termites with linear ontogeny ---------------------------------

data_So_linear = datas_So[names(datas_So) %in% c("Hsjo_So", "Kfla_So",
                                                 "Ncas_So", "PRsim_So","Znev_So")]
if (statistic) {
  dstat = data_So_linear
  data_So_linear = list(
    Hsjo_So = subset(dstat$Hsjo_So,dstat$Hsjo_So %in% linter),
    Kfla_So = subset(dstat$Kfla_So,dstat$Kfla_So %in% linter),
    Ncas_So = subset(dstat$Ncas_So,dstat$Ncas_So %in% linter),
    PRsim_So = subset(dstat$PRsim_So,dstat$PRsim_So %in% linter),
    Znev_So = subset(dstat$Znev_So,dstat$Znev_So %in% linter)
  )
  Result=supertest(data_So_linear, n = nlinter)
} else {
  Result=supertest(data_So_linear)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="So_linear"
Final_SET_Table = bind_rows(Final_SET_Table, res)
#rch = data.frame(type = "SoLin_5p",
 #                hog = return_common_hogs(res = res,degree = 5))
#Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "SoLin_4p",
                 hog = return_common_hogs(res = res,degree = 4))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)

# ------- Soldier termites with bifurcated ontogeny ---------------------------------

data_So_bifur = datas_So[names(datas_So) %in% c("Mdar_So", "Rfla_So", "Cges_So",
                                                "Mnat_So")]
if (statistic) {
  dstat = data_So_bifur
  data_So_bifur = list(
    Cges_So = subset(dstat$Cges_So,dstat$Cges_So %in% bifter),
    Mdar_So = subset(dstat$Mdar_So,dstat$Mdar_So %in% bifter),
    Mnat_So = subset(dstat$Mnat_So,dstat$Mnat_So %in% bifter),
    Rfla_So = subset(dstat$Rfla_So,dstat$Rfla_So %in% bifter)
  )
  Result=supertest(data_So_bifur, n = nbifter)
} else {
  Result=supertest(data_So_bifur)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="So_bifur"
Final_SET_Table = bind_rows(Final_SET_Table, res)
#rch = data.frame(type = "SoBif_4p",
 #                hog = return_common_hogs(res = res,degree = 4))
#Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)



# ------------------------------------------------------------------------------
# ----------- Create shared hog table ------------------------------------------
# ------------------------------------------------------------------------------

# ----------- extract table information ----------------------------------------

master$hog = as.factor(master$hog)
master = subset(master,
                subset = master$hog %in% Final_SharedHOG_table$hog,
                select = c("hog","gene","species","GO","dmel.orthologue.ID",
                           "dmel.orthologue.name","dmel.bbh.ID",
                           "amel.orthologue.ID",
                           "amel.bbh.ID",
                           "tcas.orthologue.ID",
                           "tcas.bbh.ID"))



# ----------- Merge tables -----------------------------------------------------
spelist = c()
V1 = data.frame(do.call('rbind', strsplit(as.character(gene2hog$gene), '0', fixed = T)))[,1]
V2 = data.frame(do.call('rbind', strsplit(as.character(V1), 'O', fixed = T)))[,1]
spelist = gsub("Bger", "BGER", V2)
gene2hog$spe = spelist
B = left_join(x = Final_SharedHOG_table,
              y = subset(gene2hog, spe %ni% c("Nluj")),
              by="hog")
C = left_join(x = B, y = master, by="gene")
C$hog = C$hog.x
C$hog.x = NULL
C$hog.y = NULL
D = left_join(x = C, y = SOG_list_9_termites, by="hog")


# ------------------------------------------------------------------------------
# ----------- Saving files -----------------------------------------------------
# ------------------------------------------------------------------------------

# ------------ Save the super exact test results -------------------------------

if (statistic){
  write.table(Final_SET_Table, 
              "../Result/SuperExactTest_pairwiseWOSO-stats.txt", 
              quote = F,
              row.names = F,
              sep = "\t",
              dec = ",")
}else{
  write.table(Final_SET_Table,
              "../Result/SuperExactTest_pairwiseWOSO-nostats.txt", 
              quote = F,
              row.names = F,
              sep = "\t",
              dec = ",")
}

# ------------ Save the shared hog results -------------------------------------
if (!statistic){
  write.table(D,file = "../Result/Shared_HOG_function_pairwiseWOSO_nostat.txt",
              quote = F,sep="\t",
              row.names = F,col.names = T)
}
if (statistic){
  write.table(D,file = "../Result/Shared_HOG_function_pairwiseWOSO_stat.txt",
              quote = F,sep="\t",
              row.names = F,col.names = T)
}
