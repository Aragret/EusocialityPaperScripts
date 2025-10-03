################################################################################
#
#       Title:    File preparation
#       Project:  Eusociality
#       Year:     2025
#
################################################################################


#-------------------------------------------------------------------------------
#----------- Introduction ------------------------------------------------------
#-------------------------------------------------------------------------------

'The following input files are required: 
gene2OG_map.txt
list_blue_OG.txt

'


#----------- Library installation -----------------------------------------------

library(ape)
library(MCMCglmm)
library(ggplot2)
library(cowplot)

# ---------- Set work directory ------------------------------------------------

setwd("C:/Users/caumont/Documents/WORK/Dino_McMahon/Work/Lab work and projects/10. Eusociality/MCMCglmm_CAFE/OR_count/")


# ---------- File preparation -------------------------------------------------#

OG2gene = read.table("data/gene2OG_map.txt",
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
spelist = gsub("BGER", "Bger", spelist)
spelist = gsub("KAJ", "Pame", spelist)
OG2gene$spe = spelist
OG2gene$spe = as.factor(OG2gene$spe)
summary(OG2gene)

data1 = as.data.frame.matrix(table(OG2gene$orthogroup,OG2gene$spe))
data2=data1
data2[nrow(data2) + 1,1:29] = c(colSums(data2[1:29]))
data2$Orthogroup = c(rownames(data1), "all_ors")
data2$Orthogroup = as.factor(data2$Orthogroup)
data2$Total = rowSums(x = data2[1:29])
data3 = subset(data2, data2$Total >12)
leftover13 = as.character(data3$Orthogroup)
write.table(leftover13,
            file = "data/selected_OR_OG_list.txt",
            quote=F, col.names = F, row.names = F)

# - here we added the blue branch OG --#

OG2gene = read.table("data/gene2OG_map.txt",
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
spelist = gsub("BGER", "Bger", spelist)
spelist = gsub("KAJ", "Pame", spelist)
OG2gene$spe = spelist
OG2gene$spe = as.factor(OG2gene$spe)
summary(OG2gene)

blue_OGs = read.table(file = "data/list_blue_OG.txt", header = F)
data1 = as.data.frame.matrix(table(OG2gene$orthogroup,OG2gene$spe))
data2=data1
blueOG = colSums(subset(data2, rownames(data2) %in% blue_OGs$V1))

data2[nrow(data2) + 1,1:29] = c(colSums(data2[1:29]))
data2[nrow(data2) + 1,1:29] = blueOG
data2$Orthogroup = c(rownames(data1), "all_ors", "blue_ors")
data2$Orthogroup = as.factor(data2$Orthogroup)
data2$Total = rowSums(x = data2[1:29])
data3 = subset(data2, data2$Total >12)
leftover13 = as.character(data3$Orthogroup)
write.table(leftover13,
            file = "data/selected_OR_OG_list_plusblue.txt",
            quote=F, col.names = F, row.names = F)


# ----- return to normal ----#

spe_list = c("Bger", "Dpun", "Bori", "Pame", "Cmer", "Cpun", "Mdar",
             "Hsjo", "Znev", "Kfla", "Ncas", "Cbre", "Isch", "Dlon",
             "PRsim", "Rfla", "Hten", "Cges", "Ctes", "Ssph", "Mnat",
             "Ofor", "Apac", "Aunk", "Munk", "Pred", "Ntar", "Csp4", "Nluj")


dataf = data.frame(list(Family_ID = rep("NA",dim(data3)[1]*(dim(data3)[2]-2)),
                        Species = "NA",
                        gene_count = "NA"
))

dataf$Family_ID = rep(leftover13,each = length(spe_list))
dataf$Species = rep(spe_list, dim(data3)[1])
tail(dataf,30)
spe = 'Dlon'
hog = "OG0000006"
for (hog in leftover13){
  for (spe in spe_list){
    dataf$gene_count[dataf$Family_ID == hog & dataf$Species == spe] = data3[[spe]][data3$Orthogroup==hog]
  }}



data = dataf


#------------------------ Protist -------------------------------------------#

data$Protist = "present"
spe_prot = c("Ssph", "Mnat","Ofor", "Apac", "Aunk",
             "Munk", "Pred", "Ntar", "Csp4", "Nluj")
for (spe in spe_prot) {
  data$Protist[data$Species==spe] = "absent"
}
data$Protist = factor(data$Protist)

#------------------------ Social category ------------------------------------#

data$Social_category = "NA"
soc_cat=c("Solitary",
          "Solitary",
          "Solitary",
          "Solitary",
          "Subsocial",
          "Subsocial",
          "BifurDev",
          "LinearDev",
          "LinearDev",
          "LinearDev",
          "LinearDev",
          "LinearDev",
          "LinearDev",
          "BifurDev",
          "LinearDev",
          "BifurDev",
          "BifurDev",
          "BifurDev",
          "BifurDev",
          "BifurDev",
          "BifurDev",
          "BifurDev",
          "BifurDev",
          "BifurDev",
          "BifurDev",
          "BifurDev",
          "BifurDev",
          "BifurDev",
          "BifurDev")
for (i in c(1:length(spe_list))){
  data$Social_category[data$Species==spe_list[i]]= soc_cat[i]
}

data$Social_category = factor(data$Social_category,
                              levels = c("Solitary", "Subsocial",
                                         "LinearDev", "BifurDev"))
summary(data)

# ------------------------- Food_category -----------------------------------#

data$Food_category = "NA"
food_cat=c("Litter",
           "Litter",
           "Litter",
           "Litter",
           "Wood",
           "Wood",
           "Wood",
           "Wood",
           "Wood-nesting",
           "Wood-nesting",
           "Wood-nesting",
           "Wood-nesting",
           "Wood-nesting",
           "Wood",
           "Wood",
           "Wood",
           "Wood",
           "Wood",
           "Wood",
           "Wood",
           "Fungus",
           "Fungus",
           "Humus",
           "Soil",
           "Wood",
           "Humus",
           "Humus",
           "Humus",
           "Wood")

for (i in c(1:length(spe_list))){
  data$Food_category[data$Species==spe_list[i]]= food_cat[i]
}

data$Food_category = factor(data$Food_category)
summary(data)


# ------------------------- Shortest_path ----------------------------------#

data$Shortest_path = "NA"
short_pat=c("3.000000000",
            "2.166666667",
            "3.868181818",
            "5.666666667",
            "2.333333333",
            "2.333333333",
            "2.588932806",
            "3.505494505",
            "3.571428571",
            "3.164835165",
            "3.164835165",
            "3.054945055",
            "3.472527473",
            "1.978571429",
            "2.709090909",
            "2.130434783",
            "2.071895425",
            "2.154411765",
            "2.154411765",
            "1.756876457",
            "1.717316017",
            "1.756876457",
            "1.805555556",
            "1.927272727",
            "1.424836601",
            "1.787878788",
            "2.035897436",
            "1.641025641",
            "1.845571096")

for (i in c(1:length(spe_list))){
  data$Shortest_path[data$Species==spe_list[i]]= short_pat[i]
}

data$Shortest_path = as.numeric(data$Shortest_path)
summary(data)


# ------------------------- No_terminal_caste --------------------------------#

data$No_terminal_caste = "NA"
no_tercast=c("1",
             "1",
             "1",
             "1",
             "1",
             "1",
             "4",
             "4",
             "4",
             "4",
             "4",
             "4",
             "3",
             "6",
             "3",
             "5",
             "5",
             "4",
             "4",
             "4",
             "5",
             "4",
             "2",
             "4",
             "6",
             "3",
             "4",
             "3",
             "4")

for (i in c(1:length(spe_list))){
  data$No_terminal_caste[data$Species==spe_list[i]]= no_tercast[i]
}

data$No_terminal_caste = as.numeric(data$No_terminal_caste)
summary(data)

#------------------------ Foraging -------------------------------------------#

data$Foraging = "yes"
spe_for = c("Cmer", "Cpun","Znev", "Kfla", "Ncas", "Cbre", "Isch")
for (spe in spe_for) {
  data$Foraging[data$Species==spe] = "no"
}
data$Foraging = factor(data$Foraging)
summary(data)

# ------------------------- Ontogenic_complexity_metric --------------------------#

data$Ontogenic_complexity_metric = "NA"
onto_comp=c("0.333333333",
            "0.461538462",
            "0.258519389",
            "0.176470588",
            "0.428571429",
            "0.428571429",
            "1.545038168",
            "1.141065831",
            "1.12",
            "1.263888889",
            "1.263888889",
            "1.309352518",
            "0.863924051",
            "3.032490975",
            "1.10738255",
            "2.346938776",
            "2.413249211",
            "1.85665529",
            "1.85665529",
            "2.276767945",
            "2.91152004",
            "2.276767945",
            "1.107692308",
            "2.075471698",
            "4.211009174",
            "1.677966102",
            "1.964735516",
            "1.828125",
            "2.167350805")

for (i in c(1:length(spe_list))){
  data$Ontogenic_complexity_metric[data$Species==spe_list[i]]= onto_comp[i]
}

data$Ontogenic_complexity_metric = as.numeric(data$Ontogenic_complexity_metric)
summary(data)


# ------------------------- Ontogeny -----------------------------------#

data$Ontogeny = "NA"
onto=c("linear",
       "linear",
       "linear",
       "linear",
       "linear",
       "linear",
       "bifurcated",
       "linear",
       "linear",
       "linear",
       "linear",
       "linear",
       "linear",
       "bifurcated",
       "linear",
       "bifurcated",
       "bifurcated",
       "bifurcated",
       "bifurcated",
       "bifurcated",
       "bifurcated",
       "bifurcated",
       "bifurcated",
       "bifurcated",
       "bifurcated",
       "bifurcated",
       "bifurcated",
       "bifurcated",
       "bifurcated")

for (i in c(1:length(spe_list))){
  data$Ontogeny[data$Species==spe_list[i]]= onto[i]
}

data$Ontogeny = factor(data$Ontogeny)
summary(data)


# ------------------------- Colony_size_upper_range --------------------------------#

data$Colony_size_upper_range = "NA"
col_size=c("2.4",
           "2",
           "2.6",
           "3.1",
           "1.9",
           "1.9",
           "6.8",
           "4.3",
           "4.1",
           "3.4",
           "3.8",
           "3.5",
           "3.6",
           "4",
           "3.9",
           "6.7",
           "5.5",
           "6.4",
           "6.4",
           "6",
           "6.3",
           "5.7",
           "4.6",
           "5",
           "5.2",
           "5",
           "5",
           "5",
           "5.1")

for (i in c(1:length(spe_list))){
  data$Colony_size_upper_range[data$Species==spe_list[i]]= col_size[i]
}

data$Colony_size_upper_range = as.numeric(data$Colony_size_upper_range)
summary(data)

data$gene_count = as.numeric(data$gene_count)

data$random_id = c("B_ger", "D_pun", "B_ori", "P_ame", "C_mer", "C_pun", 
                   "M_dar", "H_sjo", "Z_nev", "K_fla", "N_cas", "C_bre",
                   "I_sch", "D_lon", "PR_sim", "R_fla", "H_ten", "C_ges", 
                   "C_tes", "S_sph", "M_nat", "O_for", "A_pac", "A_unk",
                   "M_unk", "P_red", "N_tar", "C_sp4", "N_luj")


write.csv2(x = data, file = "data/dataset_mcmcglmm_ORcount.csv")
# write.csv2(x = data, file = "data/dataset_mcmcglmm_ORcount_plusblue.csv")
