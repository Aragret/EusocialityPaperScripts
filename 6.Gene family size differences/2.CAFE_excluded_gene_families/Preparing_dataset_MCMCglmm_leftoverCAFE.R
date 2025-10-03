
# This part of the script was used to make the file dataset_mcmcglmm_leftover_cafe.csv

data = read.csv("data/gene_hog_count_eusoc.csv",
                header = T,
                sep = ";",
                dec = ",", 
                stringsAsFactors = T)

str(data)
summary(data)

data$total = rowSums(data[,4:32])



HOG_pvalue = read.table("data/Gamma_all_family_results.txt", 
                        header = F)
colnames(HOG_pvalue) = c("HOG", "p-value", "signif")

leftover = as.character(data$Family_ID[!(data$Family_ID  %in% HOG_pvalue$HOG)])

data2 = subset(data, data$Family_ID %in% leftover)
data3 = subset(data2, data2$total >12)
leftover13 = as.character(data3$Family_ID)
write.table(leftover13,
            file = "data/leftover_list.txt",
            quote=F, col.names = F, row.names = F, eol = '\n')

spe_list = c("Bger", "Dpun", "Bori", "Pame", "Cmer", "Cpun", "Mdar",
             "Hsjo", "Znev", "Kfla", "Ncas", "Cbre", "Isch", "Dlon",
             "Psim", "Rfla", "Hten", "Cges", "Ctes", "Ssph", "Mnat",
             "Ofor", "Apac", "Aunk", "Munk", "Pred", "Ntar", "Cunk", "Nluj")


dataf = data.frame(list(Family_ID = rep("NA",dim(data3)[1]*(dim(data3)[2]-4)),
                        Species = "NA",
                        gene_count = "NA"
))

dataf$Family_ID = rep(leftover13,each = 29)
dataf$Species = rep(spe_list, dim(data3)[1])
tail(dataf,30)
spe = 'Dpun'
hog = "N0.HOG0000823"
for (hog in leftover13){
  for (spe in spe_list){
    dataf$gene_count[dataf$Family_ID == hog & dataf$Species == spe] = data3[[spe]][data3$Family_ID==hog]
  }}



data = dataf


#------------------------ Protist -------------------------------------------#

data$Protist = "present"
spe_prot = c("Ssph", "Mnat","Ofor", "Apac", "Aunk",
             "Munk", "Pred", "Ntar", "Cunk", "Nluj")
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

write.csv2(x = data, file = "dataset_mcmcglmm_leftover_cafe.csv")
