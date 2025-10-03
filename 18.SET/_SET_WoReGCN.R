################################################################################
#
#             Title: SET for CNG in pairwise comparison bw Wo and Re 
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

# enter the working directory.
setwd("~/WORK/Dino_McMahon/Work/Lab work and projects/6. Gene network/10.SET/3.WoReGCN/")



#-------------------------------------------------------------------------------
# ---------- Uploading annex files -------------------------------------------
#-------------------------------------------------------------------------------

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

write.table(gene2hog, "gene2hog_afterBC_inclClus2.txt",
            sep="\t",
            quote = F,
            row.names = F)

#------- Master sheet ----------------------------------------------------------

master = data.frame(read_excel("../../TABLES/small_master_sheet_eusoc.xlsx",
                               sheet = "S0_Master_sheet",
                               na = c("","NA")))


#------- list of single copy orthologues for the 9 termite species -------------

SOG_list_9_termites = data.frame(
  hog=read.table("../1.DEG/Pairwise_WORE/Script/SOG_list_9spe.txt")$V1,
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
n_linter = length(SuperExactTest::intersect(hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim))
n_bifter = length(SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Apac,hog_Cges,hog_Rfla))
nallter = length(SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Cges,hog_Rfla,
                                          hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim,hog_Apac))
allter = SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Cges,hog_Rfla,
                                   hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim,hog_Apac)
linter = SuperExactTest::intersect(hog_Hsjo,hog_Kfla,hog_Ncas,hog_PRsim)
bifter = SuperExactTest::intersect(hog_Mnat,hog_Mdar,hog_Apac,hog_Cges,hog_Rfla)


# ------------------------------------------------------------------------------
# ---------- Data-list preparation ---------------------------------------------
# --# Do not run #--------------------------------------------------------------

datas=list()

spe_list = c("Apac","Cges",
             "Hsjo","Kfla", 
             "Mdar","Mnat",
             "Ncas","PRsim",
             "Rfla",
             "Cmer","Cpun","BGER","Bori")
abg_list = c("g","b")
caste_list = c("Wo","Re")
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
  "Cryptocercus_meridianus",
  "Cryptocercus_punctulatus",
  "Blattella_germanica",
  "Blatta_orientalis")

BC = "AfterBC"
BC = "BeforeBC"

for (spe in spe_list){
  species = species_list[match(spe,spe_list)]
  file_name = paste0("data/PairwiseWoRe",BC,"_Node_groups_normalised_",
                    species,".csv")
  file_of_interest = read.csv2(file = file_name)
  exdf = data.frame(gene = file_of_interest$Node, 
                  inigroup = file_of_interest$Group_pval_Phi_tilde)
  exdf$group = file_of_interest$Group_pval_Phi_tilde
  exdf$group[exdf$inigroup == paste0("g.Network_",species,"_Er")] = paste0("g.Network_",species,"_Re")
  exdf$group[exdf$inigroup == paste0("g.Network_",species,"_Nm")] = paste0("g.Network_",species,"_Re")
  exdf$group[exdf$inigroup == paste0("g.Network_",species,"_Ad")] = paste0("g.Network_",species,"_Re")
  exdf$group[exdf$inigroup == paste0("g.Network_",species,"_Al")] = paste0("g.Network_",species,"_Re")
  exdf$group[exdf$inigroup == paste0("g.Network_",species,"_Pr")] = paste0("g.Network_",species,"_Re")
  exdf$group[exdf$inigroup == paste0("g.Network_",species,"_Ju")] = paste0("g.Network_",species,"_Wo")
  exdf$group[exdf$inigroup == paste0("b.Network_",species,"_Ju")] = paste0("b.Network_",species,"_Wo")
  if (spe %in% c("Ncas","Kfla")){
    exdf$group[exdf$inigroup == paste0("g.Network_",species,"_Al.Network_",species,"_Er")] = paste0("g.Network_",species,"_Re")
    exdf$group[exdf$inigroup == paste0("g.Network_",species,"_Al.Network_",species,"_Wo")] = "a"
    exdf$group[exdf$inigroup == paste0("g.Network_",species,"_Er.Network_",species,"_Wo")] = "a"
    exdf$group[exdf$inigroup == paste0("b.Network_",species,"_Er.Network_",species,"_Wo")] = "a"
    exdf$group[exdf$inigroup == paste0("b.Network_",species,"_Er")] = "a"
  }
  for (abg in abg_list){
    for (caste in caste_list){
      group_of_interest = paste0(abg,".Network_",species,"_", caste)
      gene_of_interest = exdf$gene[
        exdf$group == group_of_interest
      ]
      if (length(gene_of_interest) > 0){
        df = data.frame(spec= rep(spe,length(gene_of_interest)), gene=gene_of_interest)
        df1 = left_join(df, gene2hog, by= "gene")
        Temp =list(V1=unique(na.omit(df1$hog)))
        names(Temp) = paste0(spe, "_", caste,"_",abg)
        datas = append(datas, Temp)
      }
    }
  }
  print(names(datas))   

} 


# ------- Saving datas ----------------------------------------
saveRDS(datas, "Result/datas_SET_WoReGCN_AfterBC_withClus2.rds")


#saveRDS(datas, "Result/datas_SET_WoReGCN.rds")
saveRDS(datas, "Result/datas_SET_WoReGCN_BeforeBCClus2.rds")


# ------------------------------------------------------------------------------
# ------- SET and shared hog extraction ----------------------------------------
# ------------------------------------------------------------------------------

# ---------- Final table preparation -------------------------------------------

datas = readRDS("Result/datas_SET_WoReGCN_AfterBC_withClus2.rds")

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

statistic = F # F

# ------------------------------------------------------------------------------
# ------------- g --------------------------------------------------------------
# ------------------------------------------------------------------------------


# ------- All worker termites ----------------------------------------

data_Wog_termites = datas[names(datas) %in% c(
  "Apac_Wo_g","Cges_Wo_g","Hsjo_Wo_g","Kfla_Wo_g", "Mdar_Wo_g","Mnat_Wo_g",
  "Ncas_Wo_g","PRsim_Wo_g","Rfla_Wo_g")]
if (statistic) {
  dstat = data_Wog_termites
  data_Wog_termites = list(
  Apac_Wo_g = subset(dstat$Apac_Wo_g,dstat$Apac_Wo_g %in% allter),
  Cges_Wo_g = subset(dstat$Cges_Wo_g,dstat$Cges_Wo_g %in% allter),
  Hsjo_Wo_g = subset(dstat$Hsjo_Wo_g,dstat$Hsjo_Wo_g %in% allter),
  Kfla_Wo_g = subset(dstat$Kfla_Wo_g,dstat$Kfla_Wo_g %in% allter),
  Mdar_Wo_g = subset(dstat$Mdar_Wo_g,dstat$Mdar_Wo_g %in% allter),
  Mnat_Wo_g = subset(dstat$Mnat_Wo_g,dstat$Mnat_Wo_g %in% allter),
  Ncas_Wo_g = subset(dstat$Ncas_Wo_g,dstat$Ncas_Wo_g %in% allter),
  PRsim_Wo_g = subset(dstat$PRsim_Wo_g,dstat$PRsim_Wo_g %in% allter),
  Rfla_Wo_g = subset(dstat$Rfla_Wo_g,dstat$Rfla_Wo_g %in% allter)
  )
  Result=supertest(data_Wog_termites, n = nallter)
} else {
Result=supertest(data_Wog_termites)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Wog_allter"
Final_SET_Table = bind_rows(Final_SET_Table, res)

rch = data.frame(type = "Wog_9p",
                 hog = return_common_hogs(res = res,degree = 9))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Wog_7p",
                 hog = return_common_hogs(res = res,degree = 7))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Wog_5p",
                 hog = return_common_hogs(res = res,degree = 5))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Worker termites with linear ontogeny ---------------------------------
data_Wog_linear = datas[names(datas) %in% c("Hsjo_Wo_g", "Kfla_Wo_g",
                                                 "Ncas_Wo_g", "PRsim_Wo_g")]
if (statistic) {
  dstat = data_Wog_linear
  data_Wog_linear = list(
    Hsjo_Wo_g = subset(dstat$Hsjo_Wo_g,dstat$Hsjo_Wo_g %in% linter),
    Kfla_Wo_g = subset(dstat$Kfla_Wo_g,dstat$Kfla_Wo_g %in% linter),
    Ncas_Wo_g = subset(dstat$Ncas_Wo_g,dstat$Ncas_Wo_g %in% linter),
    PRsim_Wo_g = subset(dstat$PRsim_Wo_g,dstat$PRsim_Wo_g %in% linter)
  )
  Result=supertest(data_Wog_linear, n = n_linter)
} else {
  Result=supertest(data_Wog_linear)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Wog_linear"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "WoLing_4p",
                 hog = return_common_hogs(res = res,degree = 4))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Worker termites with bifurcated ontogeny ---------------------------------

data_Wog_bifur = datas[names(datas) %in% c("Mdar_Wo_g", "Rfla_Wo_g", "Cges_Wo_g",
                                                "Mnat_Wo_g", "Apac_Wo_g")]
if (statistic) {
  dstat = data_Wog_bifur
  data_Wog_bifur = list(
    Apac_Wo_g = subset(dstat$Apac_Wo_g,dstat$Apac_Wo_g %in% bifter),
    Cges_Wo_g = subset(dstat$Cges_Wo_g,dstat$Cges_Wo_g %in% bifter),
    Mdar_Wo_g = subset(dstat$Mdar_Wo_g,dstat$Mdar_Wo_g %in% bifter),
    Mnat_Wo_g = subset(dstat$Mnat_Wo_g,dstat$Mnat_Wo_g %in% bifter),
    Rfla_Wo_g = subset(dstat$Rfla_Wo_g,dstat$Rfla_Wo_g %in% bifter)
  )
  Result=supertest(data_Wog_bifur, n = n_bifter)
} else {
  Result=supertest(data_Wog_bifur)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Wog_bifur"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "WoBifg_5p",
                 hog = return_common_hogs(res = res,degree = 5))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "WoBifg_3pM",
                 hog = return_common_hogs_3M(res = res,degree = 4))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- All reproductive termites ----------------------------------------

data_Reg_termites = datas[names(datas) %in% c(
  "Cges_Re_g","Hsjo_Re_g","Kfla_Re_g", "Mdar_Re_g","Mnat_Re_g",
  "Ncas_Re_g","PRsim_Re_g","Rfla_Re_g","Apac_Re_g")]
if (statistic) {
  dstat = data_Reg_termites
  data_Reg_termites = list(
    Apac_Re_g = subset(dstat$Apac_Re_g,dstat$Apac_Re_g %in% allter),
    Cges_Re_g = subset(dstat$Cges_Re_g,dstat$Cges_Re_g %in% allter),
    Hsjo_Re_g = subset(dstat$Hsjo_Re_g,dstat$Hsjo_Re_g %in% allter),
    Kfla_Re_g = subset(dstat$Kfla_Re_g,dstat$Kfla_Re_g %in% allter),
    Mdar_Re_g = subset(dstat$Mdar_Re_g,dstat$Mdar_Re_g %in% allter),
    Mnat_Re_g = subset(dstat$Mnat_Re_g,dstat$Mnat_Re_g %in% allter),
    Ncas_Re_g = subset(dstat$Ncas_Re_g,dstat$Ncas_Re_g %in% allter),
    PRsim_Re_g = subset(dstat$PRsim_Re_g,dstat$PRsim_Re_g %in% allter),
    Rfla_Re_g = subset(dstat$Rfla_Re_g,dstat$Rfla_Re_g %in% allter)
  )
  Result=supertest(data_Reg_termites, n = nallter)
} else {
  Result=supertest(data_Reg_termites)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
# plot(Result, Layout="landscape",degree =c(8:11), sort.by="p-value", keep=FALSE,
#      show.elements=TRUE, elements.cex=0.5,cex = 0.7,
#      elements.list=subset(summary(Result)$Table,Observed.Overlap <= 20),
#      show.expected.overlap=TRUE,expected.overlap.style="box", ylim = c(0,15),
#      color.expected.overlap='blue', minMinusLog10PValue = 1.3, maxMinusLog10PValue=5,
#      color.on = rev(c("magenta")))

res$caste="Reg_allter"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "Reg_9p",
                 hog = return_common_hogs(res = res,degree = 9))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Reg_7p",
                 hog = return_common_hogs(res = res,degree = 7))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Reg_5p",
                 hog = return_common_hogs(res = res,degree = 5))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Reproductive termites with linear ontogeny ---------------------------

data_Reg_linear = datas[names(datas) %in% c("Hsjo_Re_g", "Kfla_Re_g",
                                                 "Ncas_Re_g", "PRsim_Re_g")]
if (statistic) {
  dstat = data_Reg_linear
  data_Reg_linear = list(
    Hsjo_Re_g = subset(dstat$Hsjo_Re_g,dstat$Hsjo_Re_g %in% linter),
    Kfla_Re_g = subset(dstat$Kfla_Re_g,dstat$Kfla_Re_g %in% linter),
    Ncas_Re_g = subset(dstat$Ncas_Re_g,dstat$Ncas_Re_g %in% linter),
    PRsim_Re_g = subset(dstat$PRsim_Re_g,dstat$PRsim_Re_g %in% linter)
  )
  Result=supertest(data_Reg_linear, n = n_linter)
} else {
  Result=supertest(data_Reg_linear)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Reg_linear"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "ReLing_4p",
                 hog = return_common_hogs(res = res,degree = 4))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Reproductive termites with bifurcated ontogeny -----------------------

data_Reg_bifur = datas[names(datas) %in% c("Mdar_Re_g", "Rfla_Re_g", "Cges_Re_g",
                                                "Mnat_Re_g", "Apac_Re_g")]
if (statistic) {
  dstat = data_Reg_bifur
  data_Reg_bifur = list(
    Apac_Re_g = subset(dstat$Apac_Re_g,dstat$Apac_Re_g %in% bifter),
    Cges_Re_g = subset(dstat$Cges_Re_g,dstat$Cges_Re_g %in% bifter),
    Mdar_Re_g = subset(dstat$Mdar_Re_g,dstat$Mdar_Re_g %in% bifter),
    Mnat_Re_g = subset(dstat$Mnat_Re_g,dstat$Mnat_Re_g %in% bifter),
    Rfla_Re_g = subset(dstat$Rfla_Re_g,dstat$Rfla_Re_g %in% bifter)
  )
  Result=supertest(data_Reg_bifur, n = n_bifter)
} else {
  Result=supertest(data_Reg_bifur)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Reg_bifur"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "ReBifg_5p",
                 hog = return_common_hogs(res = res,degree = 5))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "ReBifg_3pM",
                 hog = return_common_hogs_3M(res = res,degree = 4))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Juvenile subsocial cockroaches -----------------------

data_Wog_Subcockroaches = datas[names(datas) %in% c(
  "Cmer_Wo_g","Cpun_Wo_g")]
if (statistic == TRUE) {
  dstat = data_Wog_Subcockroaches
  data_Wog_Subcockroaches = list(
    Cmer_Wo_g = subset(dstat$Cmer_Wo_g,dstat$Cmer_Wo_g %in% sub_ck),
    Cpun_Wo_g = subset(dstat$Cpun_Wo_g,dstat$Cpun_Wo_g %in% sub_ck)
    )
  Result=supertest(data_Wog_Subcockroaches, n = n_sub_ck)
} else {
  Result=supertest(data_Wog_Subcockroaches)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Jug_Subck"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "JuSug_2p",
                 hog = return_common_hogs(res = res,degree = 2))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Juvenile solitary cockroaches -----------------------

data_Wog_Solcockroaches = datas[names(datas) %in% c(
  "BGER_Wo_g","Bori_Wo_g")]
if (statistic == TRUE) {
  dstat = data_Wog_Solcockroaches
  data_Wog_Solcockroaches = list(
    BGER_Wo_g = subset(dstat$BGER_Wo_g,dstat$BGER_Wo_g %in% soli_ck),
    Bori_Wo_g = subset(dstat$Bori_Wo_g,dstat$Bori_Wo_g %in% soli_ck)
  )
  Result=supertest(data_Wog_Solcockroaches, n = n_soli_ck)
} else {
  Result=supertest(data_Wog_Solcockroaches)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Jug_Solck"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "JuSog_2p",
                 hog = return_common_hogs(res = res,degree = 2))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Reproductive subsocial cockroaches -----------------------

data_Reg_Subcockroaches = datas[names(datas) %in% c(
  "Cmer_Re_g","Cpun_Re_g")]
if (statistic == TRUE) {
  dstat = data_Reg_Subcockroaches
  data_Reg_Subcockroaches = list(
    Cmer_Re_g = subset(dstat$Cmer_Re_g,dstat$Cmer_Re_g %in% sub_ck),
    Cpun_Re_g = subset(dstat$Cpun_Re_g,dstat$Cpun_Re_g %in% sub_ck)
  )
  Result=supertest(data_Reg_Subcockroaches, n = n_sub_ck)
} else {
  Result=supertest(data_Reg_Subcockroaches)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Reg_Subck"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "ReSug_2p",
                 hog = return_common_hogs(res = res,degree = 2))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Reproductive solitary cockroaches -----------------------

data_Reg_Solcockroaches = datas[names(datas) %in% c(
  "BGER_Re_g","Bori_Re_g")]
if (statistic == TRUE) {
  dstat = data_Reg_Solcockroaches
  data_Reg_Solcockroaches = list(
    BGER_Re_g = subset(dstat$BGER_Re_g,dstat$BGER_Re_g %in% soli_ck),
    Bori_Re_g = subset(dstat$Bori_Re_g,dstat$Bori_Re_g %in% soli_ck)
  )
  Result=supertest(data_Reg_Solcockroaches, n = n_soli_ck)
} else {
  Result=supertest(data_Reg_Solcockroaches)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Reg_Solck"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "ReSog_2p",
                 hog = return_common_hogs(res = res,degree = 2))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)

# ------------------------------------------------------------------------------
# ------------- b --------------------------------------------------------------
# ------------------------------------------------------------------------------

# ------- All worker termites ----------------------------------------

data_Wob_termites = datas[names(datas) %in% c(
  "Apac_Wo_b","Cges_Wo_b","Hsjo_Wo_b","Kfla_Wo_b", "Mdar_Wo_b","Mnat_Wo_b",
  "Ncas_Wo_b","PRsim_Wo_b","Rfla_Wo_b")]
if (statistic == TRUE) {
  dstat = data_Wob_termites
  data_Wob_termites = list(
    Apac_Wo_b = subset(dstat$Apac_Wo_b,dstat$Apac_Wo_b %in% allter),
    Cges_Wo_b = subset(dstat$Cges_Wo_b,dstat$Cges_Wo_b %in% allter),
    Hsjo_Wo_b = subset(dstat$Hsjo_Wo_b,dstat$Hsjo_Wo_b %in% allter),
    Kfla_Wo_b = subset(dstat$Kfla_Wo_b,dstat$Kfla_Wo_b %in% allter),
    Mdar_Wo_b = subset(dstat$Mdar_Wo_b,dstat$Mdar_Wo_b %in% allter),
    Mnat_Wo_b = subset(dstat$Mnat_Wo_b,dstat$Mnat_Wo_b %in% allter),
    Ncas_Wo_b = subset(dstat$Ncas_Wo_b,dstat$Ncas_Wo_b %in% allter),
    PRsim_Wo_b = subset(dstat$PRsim_Wo_b,dstat$PRsim_Wo_b %in% allter),
    Rfla_Wo_b = subset(dstat$Rfla_Wo_b,dstat$Rfla_Wo_b %in% allter)
  )
  Result=supertest(data_Wob_termites, n = nallter)
} else {
  Result=supertest(data_Wob_termites)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Wob_allter"
Final_SET_Table = bind_rows(Final_SET_Table, res)
#rch = data.frame(type = "Wob_9p",
 #                hog = return_common_hogs(res = res,degree = 9))
#Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
#rch = data.frame(type = "Wob_7p",
 #                hog = return_common_hogs(res = res,degree = 7))
#Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
rch = data.frame(type = "Wob_5p",
                 hog = return_common_hogs(res = res,degree = 5))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Worker termites with linear ontogeny ---------------------------------

data_Wob_linear = datas[names(datas) %in% c("Hsjo_Wo_b", "Kfla_Wo_b",
                                            "Ncas_Wo_b", "PRsim_Wo_b")]
if (statistic == TRUE) {
  dstat = data_Wob_linear
  data_Wob_linear = list(
    Hsjo_Wo_b = subset(dstat$Hsjo_Wo_b,dstat$Hsjo_Wo_b %in% linter),
    Kfla_Wo_b = subset(dstat$Kfla_Wo_b,dstat$Kfla_Wo_b %in% linter),
    Ncas_Wo_b = subset(dstat$Ncas_Wo_b,dstat$Ncas_Wo_b %in% linter),
    PRsim_Wo_b = subset(dstat$PRsim_Wo_b,dstat$PRsim_Wo_b %in% linter)
  )
  Result=supertest(data_Wob_linear, n = n_linter)
} else {
  Result=supertest(data_Wob_linear)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Wob_linear"
Final_SET_Table = bind_rows(Final_SET_Table, res)
#rch = data.frame(type = "WoLinb_4p",
 #                hog = return_common_hogs(res = res,degree = 4))
#Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Worker termites with bifurcated ontogeny ---------------------------------

data_Wob_bifur = datas[names(datas) %in% c("Mdar_Wo_b", "Rfla_Wo_b", "Cges_Wo_b",
                                           "Mnat_Wo_b", "Apac_Wo_b")]
if (statistic == TRUE) {
  dstat = data_Wob_bifur
  data_Wob_bifur = list(
    Apac_Wo_b = subset(dstat$Apac_Wo_b,dstat$Apac_Wo_b %in% bifter),
    Cges_Wo_b = subset(dstat$Cges_Wo_b,dstat$Cges_Wo_b %in% bifter),
    Mdar_Wo_b = subset(dstat$Mdar_Wo_b,dstat$Mdar_Wo_b %in% bifter),
    Mnat_Wo_b = subset(dstat$Mnat_Wo_b,dstat$Mnat_Wo_b %in% bifter),
    Rfla_Wo_b = subset(dstat$Rfla_Wo_b,dstat$Rfla_Wo_b %in% bifter)
  )
  Result=supertest(data_Wob_bifur, n = n_bifter)
} else {
  Result=supertest(data_Wob_bifur)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Wob_bifur"
Final_SET_Table = bind_rows(Final_SET_Table, res)
#rch = data.frame(type = "WoBifb_5p",
 #                hog = return_common_hogs(res = res,degree = 5))
#Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)
#rch = data.frame(type = "WoBifb_3pM",
 #                hog = return_common_hogs_3M(res = res,degree = 4))
#Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)

# ------- Juvenile subsocial cockroaches -----------------------

data_Wob_Subcockroaches = datas[names(datas) %in% c(
  "Cmer_Wo_b","Cpun_Wo_b")]
if (statistic == TRUE) {
  dstat = data_Wob_Subcockroaches
  data_Wob_Subcockroaches = list(
    Cmer_Wo_b = subset(dstat$Cmer_Wo_b,dstat$Cmer_Wo_b %in% sub_ck),
    Cpun_Wo_b = subset(dstat$Cpun_Wo_b,dstat$Cpun_Wo_b %in% sub_ck)
  )
  Result=supertest(data_Wob_Subcockroaches, n = n_sub_ck)
} else {
  Result=supertest(data_Wob_Subcockroaches)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Jub_Subck"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "JuSub_2p",
                 hog = return_common_hogs(res = res,degree = 2))
Final_SharedHOG_table = bind_rows(Final_SharedHOG_table, rch)


# ------- Juvenile solitary cockroaches -----------------------

data_Wob_Solcockroaches = datas[names(datas) %in% c(
  "BGER_Wo_b","Bori_Wo_b")]
if (statistic == TRUE) {
  dstat = data_Wob_Solcockroaches
  data_Wob_Solcockroaches = list(
    BGER_Wo_b = subset(dstat$BGER_Wo_b,dstat$BGER_Wo_b %in% soli_ck),
    Bori_Wo_b = subset(dstat$Bori_Wo_b,dstat$Bori_Wo_b %in% soli_ck)
  )
  Result=supertest(data_Wob_Solcockroaches, n = n_soli_ck)
} else {
  Result=supertest(data_Wob_Solcockroaches)}
res=data.frame(summary(Result)$Table)
if (statistic){res= res[res$Degree > 1,]}
res$caste="Jub_Solck"
Final_SET_Table = bind_rows(Final_SET_Table, res)
rch = data.frame(type = "JuSob_2p",
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
                           "dmel.bbh.name",
                           "amel.orthologue.ID",
                           "amel.bbh.ID",
                           "tcas.orthologue.ID",
                           "tcas.bbh.ID"))

master = subset(master,
                subset = master$species %in% c("Apac","Cges","Rfla",
                                               "Mnat","Mdar","Ncas",
                                               "Kfla","Hsjo","PRsim",
                                               "Cmer", "Cpun",
                                               "BGER","Bori"))


# ----------- Merge tables -----------------------------------------------------

spelist = c()
V1 = data.frame(do.call('rbind', strsplit(as.character(gene2hog$gene), '0', fixed = T)))[,1]
V2 = data.frame(do.call('rbind', strsplit(as.character(V1), 'O', fixed = T)))[,1]
spelist = gsub("Bger", "BGER", V2)
gene2hog$spe = spelist

B = left_join(x = Final_SharedHOG_table,
              y = subset(gene2hog, spe %ni% c("Znev","Nluj")),
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
              "Result/SuperExactTest_pairwise-WoRe-GCN-afterBC-stats_withClus2.txt", 
              quote = F,
              row.names = F,
              sep = "\t",
              dec = ",")
  }else{
    write.table(Final_SET_Table,
                "Result/SuperExactTest_pairwise-WoRe-GCN-afterBC-nostats_withClus2.txt", 
                quote = F,
                row.names = F,
                sep = "\t",
                dec = ",")
}

# ------------ Save the shared hog results -------------------------------------
if (!statistic){
write.table(D,file = "Result/Shared_HOG_pairwise-WoRe-GCN-AfterBC_nostat_withClus2.txt",
            quote = F,sep="\t",
            row.names = F,col.names = T)
}
if (statistic){
  write.table(D,file = "Result/Shared_HOG_pairwise-WoRe-GCN-AfterBC_stat_withClus2.txt",
              quote = F,sep="\t",
              row.names = F,col.names = T)
}

