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


# ---------- Set work directory ------------------------------------------------

setwd("C:/Users/caumont/Documents/WORK/Dino_McMahon/Work/Lab work and projects/10. Eusociality/MCMCglmm_CAFE")



print("Uploading the termite tree...")

Termite_tree <- read.tree("data/tcal_blattodea_tree.tre")

inv.phylo <- inverseA(Termite_tree, nodes = "TIPS", scale = TRUE)$Ainv

print("Uploading the termite dataset...")

data = read.csv("data/dataset_mcmcglmm_leftover_cafe.csv",
                header = T,
                sep = ";",
                dec = ",", 
                stringsAsFactors = T)


print("creating priors...")


prior.exp.phylo <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1),
            alpha.V = diag(1) * 1000)),
  R = list(V = diag(1), nu = 0.002))


prior.iw.phylo <- list(G = list(
  G1 = list(V = diag(1)*0.002/2.002, nu = 2.002)),
  R = list(V = diag(1), nu = 0.002))

#----------------------------------- Model -------------------------------------

hog="N0.HOG0000009"
hog = args[1]
print(hog)
print(paste0("Gaussian model for HOG ",hog," gene_count vs Ontogeny with expanded priors and without random IDs"))

data_subhog = subset(data, data$Family_ID == hog)
model <- MCMCglmm(gene_count ~ Ontogeny,
                  random = ~Species,
                  family = "gaussian",
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






depvar="expanded"
filenames = list.files("leftover_models", pattern = paste0("*",depvar,"*"), full.names = T) 
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
p_adj = p.adjust(p_list, method = "fdr")
df=data.frame(list(hog = hog_list,
                   p_val = p_list,
                   p_adj = p_adj))
df_psig = df[df$p_val<0.001,]
dim(df_psig)
for (no_hog in df_psig$hog){
  #no_hog = "N0.HOG0000380"
  boxplot(data[data$Family_ID==no_hog,]$gene_count~data[data$Family_ID==no_hog,]$Ontogeny, xlab = no_hog)
  points(data[data$Family_ID==no_hog,]$gene_count~data[data$Family_ID==no_hog,]$Ontogeny, col="blue",pch=19)
  points(data[data$Family_ID==no_hog & data$Species=="Mdar",]$gene_count~data[data$Family_ID==no_hog & data$Species=="Mdar",]$Ontogeny,
         pch=19, col = "red")
  points(data[data$Family_ID==no_hog & data$Species=="Dlon",]$gene_count~data[data$Family_ID==no_hog & data$Species=="Dlon",]$Ontogeny,
         pch=19, col = "green")
}
write.csv2(x = df_psig, file = "hogsize_following_ontogeny_pattern.csv",
           quote = F, row.names = F)
