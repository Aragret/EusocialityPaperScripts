################################################################################
#
#       Title:    Dataset preparation for MCMCglmm analyses
#       Project:  Eusociality
#       Year:     2025
#
################################################################################

#-------------------------------------------------------------------------------
#----------- Introduction ------------------------------------------------------
#-------------------------------------------------------------------------------

#----------- Aim ---------------------------------------------------------------

#
#The dN/dS values for Single-Copy Orthologues (SCO) and 
#Multiple-Copy Orthologues (MCO), and the relax values for SCO are uploaded.
#With this script, we add the columns for the trait values.
#

#-------------------------------------------------------------------------------
#----------- dN/dS - SCO -------------------------------------------------------
#-------------------------------------------------------------------------------


data = read.csv("data/dnds_MCO_mcmcglmm_dataset.csv",
                header = T,
                sep = ";",
                dec = ",", 
                stringsAsFactors = T)

spe_list = c("Bger", "Dpun", "Bori", "Pame", "Cmer", "Cpun", "Mdar",
             "Hsjo", "Znev", "Kfla", "Ncas", "Cbre", "Isch", "Dlon",
             "Psim", "Rfla", "Hten", "Cges", "Ctes", "Ssph", "Mnat",
             "Ofor", "Apac", "Aunk", "Munk", "Pred", "Ntar", "Cunk", "Nluj")

#------------------------ Foraging -------------------------------------------#

data$Foraging = "yes"
spe_for = c("Cmer", "Cpun","Znev", "Kfla", "Ncas", "Cbre", "Isch")
for (spe in spe_for) {
  data$Foraging[data$Species==spe] = "no"
}
data$Foraging = factor(data$Foraging)
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



write.csv2(x = data, file = "dnds_SCO.csv")








#-------------------------------------------------------------------------------
#----------- dN/dS - MCO -------------------------------------------------------
#-------------------------------------------------------------------------------



data = read.csv("data/dnds_MCO_mcmcglmm_dataset.csv",
                header = T,
                sep = ";",
                dec = ",", 
                stringsAsFactors = T)

str(data)
summary(data)

spe_list = c("Bger", "Dpun", "Bori", "Pame", "Cmer", "Cpun", "Mdar",
             "Hsjo", "Znev", "Kfla", "Ncas", "Cbre", "Isch", "Dlon",
             "Psim", "Rfla", "Hten", "Cges", "Ctes", "Ssph", "Mnat",
             "Ofor", "Apac", "Aunk", "Munk", "Pred", "Ntar", "Cunk", "Nluj")


# ------------------------ Protist --------------------------------------------#

data$Protist = "present"
spe_prot = c("Ssph", "Mnat","Ofor", "Apac", "Aunk",
             "Munk", "Pred", "Ntar", "Cunk", "Nluj")
for (spe in spe_prot) {
  data$Protist[data$Species==spe] = "absent"
}
data$Protist = factor(data$Protist)

# ------------------------ Social Categories ----------------------------------#

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







write.csv2(x = data, file = "dnds_MCO.csv")




#-------------------------------------------------------------------------------
#----------- relax - SCO -------------------------------------------------------
#-------------------------------------------------------------------------------


data = read.csv("data/relax_SCO_mcmcglmm_dataset_0.csv",
                header = T,
                sep = ";",
                dec = ",", 
                stringsAsFactors = T)

str(data)
summary(data)

spe_list = c("Bger", "Dpun", "Bori", "Pame", "Cmer", "Cpun", "Mdar",
             "Hsjo", "Znev", "Kfla", "Ncas", "Cbre", "Isch", "Dlon",
             "Psim", "Rfla", "Hten", "Cges", "Ctes", "Ssph", "Mnat",
             "Ofor", "Apac", "Aunk", "Munk", "Pred", "Ntar", "Cunk", "Nluj")


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


write.csv2(x = data, file = "relax_SCO_mcmcglmm_dataset.csv")
