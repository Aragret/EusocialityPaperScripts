#!/usr/bin/Rscript


print("importing libraries...")

library(tidyverse, lib.loc = '/home/cedria95/ca3321fu/R/x86_64-pc-linux-gnu-library/4.2')
library(ape, lib.loc = '/home/cedria95/ca3321fu/R/x86_64-pc-linux-gnu-library/4.2')
library(geiger, lib.loc = '/home/cedria95/ca3321fu/R/x86_64-pc-linux-gnu-library/4.2')
library(MCMCglmm, lib.loc = '/home/cedria95/ca3321fu/R/x86_64-pc-linux-gnu-library/4.2')
library(coda, lib.loc = '/home/cedria95/ca3321fu/R/x86_64-pc-linux-gnu-library/4.2')


print("importing arguments")
args <- commandArgs(trailingOnly = TRUE)



print("Setting working directory...")

setwd("/scratch/cedria95/Eusociality/MCMCglmm_CAFE")


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
model <- MCMCglmm(gene_count ~ Ontogeny,
                  random = ~Species,
                  family = "gaussian",
                  ginverse = list(Species = inv.phylo), 
                  prior = prior.exp.phylo,
                  data = data_subhog,
                  nitt = 13000000,
                  thin = 5000, 
                  burnin = 3000000, 
                  pr = TRUE)

print("saving...")
save(model, file = paste0("ORcount_models/Gaussian model for HOG ",hog," gene_count vs Ontogeny with expanded priors and without random IDs.Rdata"))
print("Done1")

print(paste0("Gaussian model for HOG ",hog," gene_count vs Ontogeny with iw priors and without random IDs"))

model <- MCMCglmm(gene_count ~ Ontogeny,
                  random = ~Species,
                  family = "gaussian",
                  ginverse = list(Species = inv.phylo), 
                  prior = prior.iw.phylo,
                  data = data_subhog,
                  nitt = 13000000,
                  thin = 5000, 
                  burnin = 3000000, 
                  pr = TRUE)

print("saving...")
save(model, file = paste0("ORcount_models/Gaussian model for HOG ",hog," gene_count vs Ontogeny with iw priors and without random IDs.Rdata"))
print("Done1")
