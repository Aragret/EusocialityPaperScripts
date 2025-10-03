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

'The following input files are required
gene2hog_afterBC.txt
DEG_pairwise&specif_aftBC_0.05
S0_Master_sheet
SOG_list_9spe.txt
'

#----------- Library installation ----------------------------------------------

library("RColorBrewer")
library(readxl)
library(dplyr)
library(SuperExactTest)
"%ni%" = Negate("%in%")


# ---------- Set work directory ------------------------------------------------

# enter the working directory. it must include the folder with the Node_groups
setwd("~/WORK/Dino_McMahon/Work/Lab work and projects/6. Gene network/10.SET/1.DEG/Pairwise_WORE/Script/")


# ---------- Final table preparation -------------------------------------------

datas_Wo = list()
datas_Re = list()

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

# # gene ID to HOG ID file
# gene2hog = read_xlsx(path = "../../../../TABLES/Tables_SA.xlsx",
#                      sheet = "S2_gene_to_hog",
#                      col_types = "text",
#                      skip = 2)
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

# Solitary cockroaches
hog_BGER = unique(subset(gene2hog,subset = grepl("BGER",x = gene2hog$gene))$hog)
hog_Bori = unique(subset(gene2hog,subset = grepl("Bori",x = gene2hog$gene))$hog)
n_soli_ck = length(SuperExactTest::intersect(hog_BGER,hog_Bori))
soli_ck = (SuperExactTest::intersect(hog_BGER,hog_Bori))
# subsocial cockroaches
hog_Cmer = unique(subset(gene2hog,subset = grepl("Cmer",x = gene2hog$gene))$hog)
hog_Cpun = unique(subset(gene2hog,subset = grepl("Cpun",x = gene2hog$gene))$hog)
sub_ck = (SuperExactTest::intersect(hog_Cmer,hog_Cmer))
n_sub_ck = length(SuperExactTest::intersect(hog_Cmer,hog_Cmer))
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
nlinter = length(SuperExactTest::intersect(hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim))
nbifter = length(SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Apac,hog_Cges,hog_Rfla))
nallter = length(SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Apac,hog_Cges,hog_Rfla,
                            hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim))

linter = SuperExactTest::intersect(hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim)
bifter = SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Apac,hog_Cges,hog_Rfla)
allter = SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Apac,hog_Cges,hog_Rfla,
                                            hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim)
# ---------------Initialisation ------------------------------------------------

spe_list = c("Apac","Cges","Hsjo","Kfla", "Mdar","Mnat","Ncas","PRsim","Rfla",
             "Znev","Cmer","Cpun","Bger","Bori")

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

colu_wore = c("Pr_Wo","Al_Wo","Ne_Wo","Er_Wo","Nm_Wo","Ad_Ju", "Re_Wo")

# ------------------------------------------------------------------------------

for (spe in spe_list){
  species = species_list[match(spe,spe_list)]
  print(spe)
  data_sub = subset(data_raw, subset= data_raw$Species == spe)
  sub_sub_data = subset(data_sub, select = colnames(data_sub) %in% c("gene_name",colu_wore))
  sub_sub_data$Re_Wo = NA
  if (spe %in% c("Ncas","Kfla")){
    NK = sub_sub_data
    NK$Re_Wo[NK$Al_Wo == "none" & NK$Ne_Wo == "none"] ="none"
    NK$Re_Wo[NK$Al_Wo == "Al" & NK$Ne_Wo == "Ne"] = "Re"
    NK$Re_Wo[NK$Al_Wo == "Wo" & NK$Ne_Wo == "Wo"] = "Wo"
    NK$Re_Wo[NK$Al_Wo == "Al" & NK$Ne_Wo == "Wo"] = "none"
    NK$Re_Wo[NK$Al_Wo == "Wo" & NK$Ne_Wo == "Ne"] = "none"
    NK$Re_Wo[NK$Al_Wo == "Al" & NK$Ne_Wo == "none"] = "Re"
    NK$Re_Wo[NK$Al_Wo == "none" & NK$Ne_Wo == "Ne"] ="Re"
    NK$Re_Wo[NK$Al_Wo == "Wo" & NK$Ne_Wo == "none"] ="Wo"
    NK$Re_Wo[NK$Al_Wo == "none" & NK$Ne_Wo == "Wo"]="Wo"
    sub_sub_data = NK
  }
  for (colu in colu_wore){
    if ((spe %ni% c("Ncas","Kfla") & !is.na(sub_sub_data[1,colu]))| (spe %in% c("Ncas","Kfla") & colu== "Re_Wo")){
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
        if (casteA %in% c("Wo","Ju")){
          names(Temp) = paste0(spe,"_Wo")
          datas_Wo = append(datas_Wo, Temp)
        }
        if (casteA %ni% c("Wo","Ju")){
          names(Temp) = paste0(spe,"_Re")
          datas_Re = append(datas_Re, Temp)
        }
      }
    }       
  }
}

saveRDS(datas_Re, "../Result/data_biased_Re_aftBC.rds")
saveRDS(datas_Wo, "../Result/data_biased_Wo_aftBC.rds")
names(datas_Wo)

# ------------------------------------------------------------------------------
# ------- SET and shared hog extraction ----------------------------------------
# ------------------------------------------------------------------------------
datas_Wo = readRDS("../Result/data_biased_Wo_aftBC.rds")
datas_Re = readRDS("../Result/data_biased_Re_aftBC.rds")
statistic = F

# ------- All reproductive termites ----------------------------------------

data_Re_termites = datas_Re[names(datas_Re) %ni% c("Cmer_Re", "Cpun_Re",
                                                   "Bger_Re", "Bori_Re")]
names(data_Re_termites)
if (statistic) {
  dstat = data_Re_termites
  data_Re_termites = list(
    Apac_Re = subset(dstat$Apac_Re,dstat$Apac_Re %in% allter),
    Cges_Re = subset(dstat$Cges_Re,dstat$Cges_Re %in% allter),
    Hsjo_Re = subset(dstat$Hsjo_Re,dstat$Hsjo_Re %in% allter),
    Kfla_Re = subset(dstat$Kfla_Re,dstat$Kfla_Re %in% allter),
    Mdar_Re = subset(dstat$Mdar_Re,dstat$Mdar_Re %in% allter),
    Mnat_Re = subset(dstat$Mnat_Re,dstat$Mnat_Re %in% allter),
    Ncas_Re = subset(dstat$Ncas_Re,dstat$Ncas_Re %in% allter),
    PRsim_Re = subset(dstat$PRsim_Re,dstat$PRsim_Re %in% allter),
    Rfla_Re = subset(dstat$Rfla_Re,dstat$Rfla_Re %in% allter)
  )
  Result=supertest(data_Re_termites, n = nallter)
} else {
  Result=supertest(data_Re_termites)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Re_allter"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "Re_9p",
                hog = return_common_hogs(res = res,degree = 9))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Re_7p",
                hog = return_common_hogs(res = res,degree = 7))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Re_5p",
                hog = return_common_hogs(res = res,degree = 5))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- All Worker termites ----------------------------------------

data_Wo_termites = datas_Wo[names(datas_Wo) %ni% c("Cmer_Wo", "Cpun_Wo", 
                                                   "Bger_Wo", "Bori_Wo")]
if (statistic) {
  dstat = data_Wo_termites
  data_Wo_termites = list(
    Apac_Wo = subset(dstat$Apac_Wo,dstat$Apac_Wo %in% allter),
    Cges_Wo = subset(dstat$Cges_Wo,dstat$Cges_Wo %in% allter),
    Hsjo_Wo = subset(dstat$Hsjo_Wo,dstat$Hsjo_Wo %in% allter),
    Kfla_Wo = subset(dstat$Kfla_Wo,dstat$Kfla_Wo %in% allter),
    Mdar_Wo = subset(dstat$Mdar_Wo,dstat$Mdar_Wo %in% allter),
    Mnat_Wo = subset(dstat$Mnat_Wo,dstat$Mnat_Wo %in% allter),
    Ncas_Wo = subset(dstat$Ncas_Wo,dstat$Ncas_Wo %in% allter),
    PRsim_Wo = subset(dstat$PRsim_Wo,dstat$PRsim_Wo %in% allter),
    Rfla_Wo = subset(dstat$Rfla_Wo,dstat$Rfla_Wo %in% allter)
  )
  Result=supertest(data_Wo_termites, n = nallter)
} else {
  Result=supertest(data_Wo_termites)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Wo_allter"
Final_SET_Table = bind_rows(Final_SET_Table, res)
#rch = data.frame(type = "Wo_9p",
 #              hog = return_common_hogs(res = res,degree = 9)) 
#Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Wo_7p",
                hog = return_common_hogs(res = res,degree = 7))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Wo_5p",
                hog = return_common_hogs(res = res,degree = 5))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Reproductive termites with linear ontogeny ---------------------------

data_Re_linear = datas_Re[names(datas_Re) %in% c("Hsjo_Re", "Kfla_Re",
                                                 "Ncas_Re", "PRsim_Re")]
if (statistic) {
  dstat = data_Re_linear
  data_Re_linear = list(
    Hsjo_Re = subset(dstat$Hsjo_Re,dstat$Hsjo_Re %in% linter),
    Kfla_Re = subset(dstat$Kfla_Re,dstat$Kfla_Re %in% linter),
    Ncas_Re = subset(dstat$Ncas_Re,dstat$Ncas_Re %in% linter),
    PRsim_Re = subset(dstat$PRsim_Re,dstat$PRsim_Re %in% linter)
  )
  Result=supertest(data_Re_linear, n = nlinter)
} else {
  Result=supertest(data_Re_linear)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Re_linear"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "ReLin_4p",
                hog = return_common_hogs(res = res,degree = 4))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Reproductive termites with bifurcated ontogeny -----------------------

data_Re_bifur = datas_Re[names(datas_Re) %in% c("Mdar_Re", "Rfla_Re", "Cges_Re",
                                                "Mnat_Re", "Apac_Re")]
if (statistic) {
  dstat = data_Re_bifur
  data_Re_bifur = list(
    Apac_Re = subset(dstat$Apac_Re,dstat$Apac_Re %in% bifter),
    Cges_Re = subset(dstat$Cges_Re,dstat$Cges_Re %in% bifter),
    Mdar_Re = subset(dstat$Mdar_Re,dstat$Mdar_Re %in% bifter),
    Mnat_Re = subset(dstat$Mnat_Re,dstat$Mnat_Re %in% bifter),
    Rfla_Re = subset(dstat$Rfla_Re,dstat$Rfla_Re %in% bifter)
  )
  Result=supertest(data_Re_bifur, n = nbifter)
} else {
  Result=supertest(data_Re_bifur)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Re_bifur"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "ReBif_5p",
                hog = return_common_hogs(res = res,degree = 5))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "ReBif_3pM",
                hog = return_common_hogs_3M(res = res,degree = 4))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Worker termites with linear ontogeny ---------------------------------

data_Wo_linear = datas_Wo[names(datas_Wo) %in% c("Hsjo_Wo", "Kfla_Wo",
                                                 "Ncas_Wo", "PRsim_Wo")]
if (statistic) {
  dstat = data_Wo_linear
  data_Wo_linear = list(
    Hsjo_Wo = subset(dstat$Hsjo_Wo,dstat$Hsjo_Wo %in% linter),
    Kfla_Wo = subset(dstat$Kfla_Wo,dstat$Kfla_Wo %in% linter),
    Ncas_Wo = subset(dstat$Ncas_Wo,dstat$Ncas_Wo %in% linter),
    PRsim_Wo = subset(dstat$PRsim_Wo,dstat$PRsim_Wo %in% linter)
  )
  Result=supertest(data_Wo_linear, n = nlinter)
} else {
  Result=supertest(data_Wo_linear)}
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Wo_linear"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "WoLin_4p",
                hog = return_common_hogs(res = res,degree = 4))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Worker termites with bifurcated ontogeny ---------------------------------

data_Wo_bifur = datas_Wo[names(datas_Wo) %in% c("Mdar_Wo", "Rfla_Wo", "Cges_Wo",
                                                "Mnat_Wo", "Apac_Wo")]
if (statistic) {
  dstat = data_Wo_bifur
  data_Wo_bifur = list(
    Apac_Wo = subset(dstat$Apac_Wo,dstat$Apac_Wo %in% bifter),
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
rch = data.frame(type = "WoBif_5p",
                hog = return_common_hogs(res = res,degree = 5))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "WoBif_3pM",
                hog = return_common_hogs_3M(res = res,degree = 4))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- CK JUSU ---------------------------------

data_Wo_JUSU = datas_Wo[names(datas_Wo) %in% c("Cpun_Wo", "Cmer_Wo")]
dstat = data_Wo_JUSU
data_Wo_JUSU = list(
    Cpun_Wo = subset(dstat$Cpun_Wo,dstat$Cpun_Wo %in% sub_ck),
    Cmer_Wo = subset(dstat$Cmer_Wo,dstat$Cmer_Wo %in% sub_ck)
)
Result=supertest(data_Wo_JUSU, n = n_sub_ck)
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Wo_JUSU"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "Wo_JUSU",
                 hog = return_common_hogs(res = res,degree = 2))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- CK JUSO ---------------------------------

data_Wo_JUSO = datas_Wo[names(datas_Wo) %in% c("Bger_Wo", "Bori_Wo")]
dstat = data_Wo_JUSO
data_Wo_JUSO = list(
  Bger_Wo = subset(dstat$Bger_Wo,dstat$Bger_Wo %in% soli_ck),
  Bori_Wo = subset(dstat$Bori_Wo,dstat$Bori_Wo %in% soli_ck)
)
Result=supertest(data_Wo_JUSO, n = n_soli_ck)
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Wo_JUSO"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "Wo_JUSO",
                 hog = return_common_hogs(res = res,degree = 2))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- CK ADSU ---------------------------------

data_Wo_bifur = datas_Re[names(datas_Re) %in% c("Cmer_Re", "Cpun_Re")]
dstat = data_Wo_bifur
data_Wo_bifur = list(
  Cpun_Re = subset(dstat$Cpun_Re,dstat$Cpun_Re %in% sub_ck),
  Cmer_Re = subset(dstat$Cmer_Re,dstat$Cmer_Re %in% sub_ck)
)
Result=supertest(data_Wo_bifur, n = n_sub_ck)
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Ad_ADSU"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "Ad_ADSU",
                 hog = return_common_hogs(res = res,degree = 2))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)

# ------- CK ADSO ---------------------------------

data_Wo_bifur = datas_Re[names(datas_Re) %in% c("Bger_Re", "Bori_Re")]
dstat = data_Wo_bifur
data_Wo_bifur = list(
  Bger_Re = subset(dstat$Bger_Re,dstat$Bger_Re %in% soli_ck),
  Bori_Re = subset(dstat$Bori_Re,dstat$Bori_Re %in% soli_ck)
)
Result=supertest(data_Wo_bifur, n = n_soli_ck)
res=data.frame(summary(Result)$Table)
res= res[res$Degree > 1,]
res$caste="Ad_ADSO"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "Ad_ADSO",
                 hog = return_common_hogs(res = res,degree = 2))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


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
              "../Result/SuperExactTest_pairwiseWORE-stats.txt", 
              quote = F,
              row.names = F,
              sep = "\t",
              dec = ",")
}else{
  write.table(Final_SET_Table,
              "../Result/SuperExactTest_pairwiseWORE-nostats.txt", 
              quote = F,
              row.names = F,
              sep = "\t",
              dec = ",")
}

# ------------ Save the shared hog results -------------------------------------
if (!statistic){
  write.table(D,file = "../Result/Shared_HOG_function_pairwiseWORE_nostat.txt",
              quote = F,sep="\t",
              row.names = F,col.names = T)
}
if (statistic){
  write.table(D,file = "../Result/Shared_HOG_function_pairwiseWORE_stat.txt",
              quote = F,sep="\t",
              row.names = F,col.names = T)}




# ------------ Shared with Cockroaches ---------------------------------------

datas_Wo = readRDS("../Result/data_biased_Wo_aftBC.rds")
datas_Re = readRDS("../Result/data_biased_Re_aftBC.rds")


scockroaches = data.frame(read_excel("../Result/Results_pDEG_WORE.xlsx",
                               sheet = "Shared with Cockroaches",
                               na = c("","NA"),skip = 1))


scockroaches$Bger.Ju[scockroaches$hog %in% datas_Wo$Bger_Wo] = "Yes"
scockroaches$Bori.Ju[scockroaches$hog %in% datas_Wo$Bori_Wo] = "Yes"
scockroaches$Cmer.Ju[scockroaches$hog %in% datas_Wo$Cmer_Wo] = "Yes"
scockroaches$Cpun.Ju[scockroaches$hog %in% datas_Wo$Cpun_Wo] = "Yes"
scockroaches$Bger.Ad[scockroaches$hog %in% datas_Re$Bger_Re] = "Yes"
scockroaches$Bori.Ad[scockroaches$hog %in% datas_Re$Bori_Re] = "Yes"
scockroaches$Cmer.Ad[scockroaches$hog %in% datas_Re$Cmer_Re] = "Yes"
scockroaches$Cpun.Ad[scockroaches$hog %in% datas_Re$Cpun_Re] = "Yes"

write.csv2(scockroaches, "../Result/Shared with Cockroaches - pDEG_WORE.csv",
           row.names = F)
