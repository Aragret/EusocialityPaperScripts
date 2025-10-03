################################################################################
#
#       Title:    MCMCglmm_leftover_CAFE
#       Project:  Eusociality
#       Year:     2025
#
################################################################################


#-------------------------------------------------------------------------------
#----------- Introduction ------------------------------------------------------
#-------------------------------------------------------------------------------

#----------- Library installation -----------------------------------------------

library(ape)
library(MCMCglmm)
library(ggplot2)
library(cowplot)

# ---------- Set work directory ------------------------------------------------

setwd("C:/Users/caumont/Documents/WORK/Dino_McMahon/Work/Lab work and projects/10. Eusociality/MCMCglmm_CAFE/OR_count/")

##### do not run ####
# ---------- File preparation -------------------------------------------------#
# This part of the script was used to make the file dataset_mcmcglmm_leftover_cafe.csv

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

#####
# Ran on HPC#
print("Uploading the termite tree...")

Termite_tree <- read.tree("data/tcal_blattodea_tree.tre")
Termite_tree$tip.label <- gsub("Psim", "PRsim", Termite_tree$tip.label)
Termite_tree$tip.label <- gsub("Cunk", "Csp4", Termite_tree$tip.label)

inv.phylo <- inverseA(Termite_tree, nodes = "TIPS", scale = TRUE)$Ainv

print("Uploading the termite dataset...")

data = read.csv("data/dataset_mcmcglmm_ORcount.csv",
                header = T,
                sep = ";",
                dec = ",", 
                stringsAsFactors = T)

str(data)
print("creating priors...")


prior.exp.phylo <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1),
            alpha.V = diag(1) * 1000)),
  R = list(V = diag(1), nu = 0.002))


prior.iw.phylo <- list(G = list(
  G1 = list(V = diag(1)*0.002/2.002, nu = 2.002)),
  R = list(V = diag(1), nu = 0.002))


#----------------------------------- Model -------------------------------------


hog = args[1]
print(hog)
print(paste0("Gaussian model for HOG ",hog," gene_count vs Ontogeny with expanded priors and without random IDs"))

data_subhog = subset(data, data$Family_ID == hog)
data_subhog = subset(x = data_subhog,subset = Social_category %in% c("BifurDev", "LinearDev"))
modelalle <- MCMCglmm(gene_count ~ Ontogeny,
                  random = ~Species,
                  family = "exponential",
                  ginverse = list(Species = inv.phylo), 
                  prior = prior.exp.phylo,
                  data = data_subhog,
                  nitt = 1300000,
                  thin = 500, 
                  burnin = 300000, 
                  pr = TRUE)

print("saving...")
save(model, file = paste0("Gaussian model for HOG ",hog," gene_count vs Ontogeny with expanded priors and without random IDs.Rdata"))
print("Done1")

print(paste0("Gaussian model for HOG ",hog," gene_count vs Ontogeny with iw priors and without random IDs"))

model <- MCMCglmm(gene_count ~ Ontogeny,
                  random = ~Species,
                  family = "gaussian",
                  ginverse = list(Species = inv.phylo), 
                  prior = prior.iw.phylo,
                  data = data_subhog,
                  nitt = 1300000,
                  thin = 500, 
                  burnin = 300000, 
                  pr = TRUE)

print("saving...")
save(model, file = paste0("Gaussian model for HOG ",hog," gene_count vs Ontogeny with iw priors and without random IDs.Rdata"))
print("Done1")


#END Ran on HPC END#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ---------- The Total OR model ------------------------------------------------


load("ORcount_models_totals/Gaussian model for HOG all_ors gene_count vs Ontogeny with expanded priors and with random IDs.Rdata")
summary(model)$solutions
#                post.mean  l-95% CI    u-95% CI eff.samp  pMCMC
# (Intercept)    110.80928  78.27488 137.8080805 1850.405 0.0005
# Ontogenylinear -23.11214 -44.58965  -0.7210318 2000.000 0.0390
load("ORcount_models_totals/Gaussian model for HOG all_ors gene_count vs Social_category with expanded priors and with random IDs.Rdata")
summary(model)$solutions
#                           post.mean  l-95% CI   u-95% CI eff.samp  pMCMC
# (Intercept)              101.368139  66.82976 130.473672     2000 0.0005
# Social_categoryLinearDev -29.656087 -49.33358  -8.458558     2000 0.0080
# Social_categorySolitary   -5.861714 -43.53376  34.529795     2000 0.7380
# Social_categorySubsocial  -1.771242 -39.82725  38.127836     2000 0.9270

load("ORcount_models_OCM/Gaussian model for HOG all_ors gene_count vs Ontogenic_complexity_metric with expanded priors and with random IDs.Rdata")
summary(model)$solutions
#                             post.mean  l-95% CI  u-95% CI eff.samp  pMCMC
# (Intercept)                 91.656962 62.174302 116.39733     2000 0.0005
# Ontogenic_complexity_metric  0.849257 -8.385558  11.01013     2000 0.8690
# All tests gave the same answer: (except exp without model not significant for Lin vs Bif species) 
# significant difference between Linear and bif and between tbif and tlin
# no significance with OCM

#Figure
library(ggplot2)
library(cowplot)
pdf("Result/total gene count  - quick figure updated.pdf")

no_hog = "all_ors"
ggdata = data[data$Family_ID==no_hog,]
ggdata$colour = "Others"
ggdata$colour[ggdata$Species == "Mdar"] = "Mdar"
ggdata$colour[ggdata$Species == "Dlon"] = "Dlon"
Bin = ggplot(data =ggdata )+
  geom_jitter(aes(Ontogeny, gene_count, color = colour),width = .1,height = 0)+
  scale_color_manual(name = "Species", 
                     values = c("Others" = "grey", "Mdar" = "blue","Dlon" = "green"))+
  theme_bw()+
  ggtitle(no_hog)

Soc = ggplot(data =ggdata )+
  geom_jitter(aes(Social_category, gene_count, color = colour),width = .1,height = 0)+
  scale_color_manual(name = "Species", 
                     values = c("Others" = "grey", "Mdar" = "blue","Dlon" = "green"))+
  theme_bw()+
  ggtitle(no_hog)  

plot_grid(Bin, Soc, labels = "AUTO")

ggsave2(filename = "Result/total gene count  - quick figure updated.pdf",
        width = 10, height = 6)

# ------------ checking the OrHOGs ---------------------------------------------


fin_df = data.frame(hog = character(),
                    p_val = numeric(),
                    p_adj = numeric(),
                    covariate = character(),
                    prior = character(),
                    randon = character())
for (var1 in c("Ontogeny", "Social_category")){
for (prior in c("iw", "expanded")){
for (rd in c("with", "without")){
depvar=paste0(var1," with ",prior," priors and ",rd," random")
filenames = list.files("ORcount_models", pattern = paste0("*",depvar,"*"), full.names = T) 
N_row = length(filenames)
j=0
p_list=c()
hog_list=c()
for (modfile in filenames){
  print(modfile)
  j = j + 1
  mod = load(modfile)
  model = eval(as.name(mod))
  hog_name = strsplit(modfile, " ")[[1]][5]
  hog_list=append(hog_list, hog_name)
  p_temp = summary(model)$solutions[2,5]
  p_list=append(p_list, p_temp)
}
p_adj = p.adjust(p_list,method = "fdr")
df=data.frame(list(hog = hog_list,
                   p_val = p_list,
                   p_adj = p_adj,
                   covariate = var1,
                   prior = prior,
                   randon = rd))
df_psig = df[df$p_adj<0.05,]
dim(df_psig)
fin_df = bind_rows(fin_df,df_psig)
}}}

write.csv2(x = fin_df, file = paste0("Result/Significance OrHOG_possvm01.csv"),
           quote = F, row.names = F)


# ----------- Quick figures ----------------------------------------------------


# Ontogeny

depvar="Ontogeny with iw priors and with random"
filenames = list.files("ORcount_models", pattern = paste0("*",depvar,"*"), full.names = T) 
N_row = length(filenames)
j=0
p_list=c()
hog_list=c()
for (modfile in filenames){
  print(modfile)
  j = j + 1
  mod = load(modfile)
  model = eval(as.name(mod))
  hog_name = strsplit(modfile, " ")[[1]][5]
  hog_list=append(hog_list, hog_name)
  p_temp = summary(model)$solutions[2,5]
  p_list=append(p_list, p_temp)
}
p_adj = p.adjust(p_list,method = "fdr")
df=data.frame(list(hog = hog_list,
                   p_val = p_list,
                   p_adj = p_adj))
df_psig = df[df$p_adj<0.05,]
dim(df_psig)

pdf("Result/sigOrHOG ontogeny  - quick figure.pdf",width = 6,height = 9)
par(mfrow=c(3,3))
for (no_hog in df_psig$hog){
  #no_hog = "N0.HOG0000380"
  boxplot(data[data$Family_ID==no_hog,]$gene_count~data[data$Family_ID==no_hog,]$Ontogeny, xlab = no_hog)
  points(data[data$Family_ID==no_hog,]$gene_count~data[data$Family_ID==no_hog,]$Ontogeny, col="blue",pch=19)
  points(data[data$Family_ID==no_hog & data$Species=="Mdar",]$gene_count~data[data$Family_ID==no_hog & data$Species=="Mdar",]$Ontogeny,
         pch=19, col = "red")
  points(data[data$Family_ID==no_hog & data$Species=="Dlon",]$gene_count~data[data$Family_ID==no_hog & data$Species=="Dlon",]$Ontogeny,
         pch=19, col = "green")
}
dev.off()


# Social_category

depvar="Social_category with iw priors and with random"
filenames = list.files("ORcount_models", pattern = paste0("*",depvar,"*"), full.names = T) 
N_row = length(filenames)
j=0
p_list=c()
hog_list=c()
for (modfile in filenames){
  print(modfile)
  j = j + 1
  mod = load(modfile)
  model = eval(as.name(mod))
  hog_name = strsplit(modfile, " ")[[1]][5]
  hog_list=append(hog_list, hog_name)
  p_temp = summary(model)$solutions[2,5]
  p_list=append(p_list, p_temp)
}
p_adj = p.adjust(p_list,method = "fdr")
df=data.frame(list(hog = hog_list,
                   p_val = p_list,
                   p_adj = p_adj))
df_psig = df[df$p_adj<0.05,]
dim(df_psig)


pdf("Result/sigOrHOG Social_category  - quick figure.pdf",
    width = 10,height = 6)
par(mfrow=c(2,3))
for (no_hog in unique(fin_df$hog)){
  #no_hog = "OG119"
  boxplot(data[data$Family_ID==no_hog,]$gene_count~data[data$Family_ID==no_hog,]$Social_category, xlab = no_hog)
  points(data[data$Family_ID==no_hog,]$gene_count~data[data$Family_ID==no_hog,]$Social_category, col="blue",pch=19)
  points(data[data$Family_ID==no_hog & data$Species=="Mdar",]$gene_count~data[data$Family_ID==no_hog & data$Species=="Mdar",]$Social_category,
         pch=19, col = "red")
  points(data[data$Family_ID==no_hog & data$Species=="Dlon",]$gene_count~data[data$Family_ID==no_hog & data$Species=="Dlon",]$Social_category,
         pch=19, col = "green")
}
for (no_hog in unique(fin_df$hog)){
  #no_hog = "OG119"
  boxplot(data[data$Family_ID==no_hog,]$gene_count~data[data$Family_ID==no_hog,]$Ontogeny, xlab = no_hog)
  points(data[data$Family_ID==no_hog,]$gene_count~data[data$Family_ID==no_hog,]$Ontogeny, col="blue",pch=19)
  points(data[data$Family_ID==no_hog & data$Species=="Mdar",]$gene_count~data[data$Family_ID==no_hog & data$Species=="Mdar",]$Ontogeny,
         pch=19, col = "red")
  points(data[data$Family_ID==no_hog & data$Species=="Dlon",]$gene_count~data[data$Family_ID==no_hog & data$Species=="Dlon",]$Ontogeny,
         pch=19, col = "green")
}
dev.off()



# info for main text

datablue = read.csv("data/dataset_mcmcglmm_ORcount_plusblue.csv",
                    header = T,
                    sep = ";",
                    dec = ",", 
                    stringsAsFactors = T)

# Mean number of ORs (cluster 1) per TW species  
mean(datablue[datablue$Family_ID=="OG112" & datablue$Social_category == "BifurDev",]$gene_count)
# 20.9375

# Mean number of ORs (cluster 1) per FW species 
mean(datablue[datablue$Family_ID=="OG112" & datablue$Social_category == "LinearDev",]$gene_count)
# 10.42857

# MCMCglmm output Cluster 1
load("ORcount_models/Gaussian model for HOG OG112 gene_count vs Social_category with expanded priors and with random IDs.Rdata")
summary(model)$solutions
#                          post.mean  l-95% CI  u-95% CI eff.samp pMCMC
# (Intercept)               23.19384  12.07170 35.697959 2000.000 0.001
# Social_categoryLinearDev -10.08762 -17.21458 -2.551893 1778.064 0.007
# Social_categorySolitary  -22.50572 -36.15525 -7.954913 2000.000 0.005
# Social_categorySubsocial -10.72234 -25.09009  3.738953 1813.485 0.142


# Mean number of ORs (cluster 2) per TW species  
mean(datablue[datablue$Family_ID=="blue_ors" & datablue$Social_category == "BifurDev",]$gene_count)
# 25.875

# Mean number of ORs (cluster 2) per FW species 
mean(datablue[datablue$Family_ID=="blue_ors" & datablue$Social_category == "LinearDev",]$gene_count)
# 3


# MCMCglmm output Cluster 2
load("ORcount_models/Gaussian model for HOG blue_ors gene_count vs Social_category with expanded priors and with random IDs.Rdata")
summary(model)$solutions
#                          post.mean   l-95% CI   u-95% CI eff.samp pMCMC
# (Intercept)               20.27527   4.143555 33.5838117 2000.000 0.026
# Social_categoryLinearDev -17.90136 -28.930707 -6.9446062 2000.000 0.008
# Social_categorySolitary  -20.31562 -38.659266 -0.0663502 1511.727 0.054
# Social_categorySubsocial -13.14948 -31.659413  8.0318592 2000.000 0.201


# Mean number of total ORs per TW species  
mean(datablue[datablue$Family_ID=="all_ors" & datablue$Social_category == "BifurDev",]$gene_count)
# 104.6875

# Mean number of total ORs per FW species 
mean(datablue[datablue$Family_ID=="all_ors" & datablue$Social_category == "LinearDev",]$gene_count)
# 72.71429


load("ORcount_models_totals/Gaussian model for HOG all_ors gene_count vs Social_category with expanded priors and with random IDs.Rdata")
summary(model)$solutions
#                           post.mean  l-95% CI   u-95% CI eff.samp  pMCMC
# (Intercept)              101.368139  66.82976 130.473672     2000 0.0005
# Social_categoryLinearDev -29.656087 -49.33358  -8.458558     2000 0.0080
# Social_categorySolitary   -5.861714 -43.53376  34.529795     2000 0.7380
# Social_categorySubsocial  -1.771242 -39.82725  38.127836     2000 0.9270




# Figure 

# no_hog = "OG112"
#no_hog = "OG202"
#no_hog = "OG21"
# no_hog = "blue_ors"
no_hog = "all_ors"
ggdata = datablue[datablue$Family_ID==no_hog,]
ggdata$colour[ggdata$Social_category == "Solitary"] = "Solitary"
ggdata$colour[ggdata$Social_category == "Subsocial"] = "Subsocial"
ggdata$colour[ggdata$Social_category == "LinearDev"] = "Social (FW)"
ggdata$colour[ggdata$Social_category == "BifurDev"] = "Social (TW)"
ggdata$colour[ggdata$Species == "Mdar"] = "Mdar"
ggdata$colour[ggdata$Species == "Dlon"] = "Dlon"
ggdata$Social_category = factor(ggdata$Social_category,
                                levels = c("Solitary", "Subsocial", 
                                "LinearDev" , "BifurDev"),
                                labels =  c("Solitary", "Subsocial",
                                            "Social (FW)", "Social (TW)")
                                )
ggdata$Ontogeny = factor(ggdata$Ontogeny, levels = c("linear", "bifurcated"))
Bin = ggplot(data =ggdata )+
  geom_jitter(aes(Ontogeny, gene_count, color = colour),width = .1,height = 0)+
  scale_color_manual(name = "Species", 
                     values = c("Solitary" = "grey", 
                                "Subsocial" = "grey50",
                                "Social (FW)" = "black",
                                "Social (TW)" = "deepskyblue",
                                "Mdar" = "red","Dlon" = "green"))+
  theme_half_open()+
  ggtitle(no_hog)+ ylab("Gene count")
Soc = ggplot(data =ggdata )+
  geom_jitter(aes(Social_category, gene_count, color = colour),width = .2,height = 0)+
  scale_color_manual(name = "Species", 
                     values = c("Solitary" = "grey", 
                                "Subsocial" = "grey50",
                                "Social (FW)" = "black",
                                "Social (TW)" = "deepskyblue",
                                "Mdar" = "red","Dlon" = "green"))+
  theme_half_open()+
  ggtitle(no_hog)+ ylab("Gene count")
plot_grid(Bin, Soc, labels = "AUTO")

ggsave2(filename = paste0("Result/",no_hog," gene count - quick figure 28-07.pdf"),
        width = 14, height = 6)

