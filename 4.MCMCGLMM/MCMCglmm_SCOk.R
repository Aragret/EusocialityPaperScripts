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

setwd("/scratch/cedria95/Eusociality/MCMCglmm")


print("Uploading the termite tree...")

Termite_tree <- read.tree("data/tcal_blattodea_tree.tre")

inv.phylo <- inverseA(Termite_tree, nodes = "TIPS", scale = TRUE)$Ainv

print("Uploading the termite dataset...")

data = read.csv("data/relax_SCO_mcmcglmm_dataset.csv",
                header = T,
                sep = ";",
                dec = ",", 
                stringsAsFactors = T)

names(data)
str(data)

print("creating priors...")


prior.exp.phylo <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1),
            alpha.V = diag(1) * 1000),
  G2 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), 
            alpha.V = diag(1) * 1000)),
  R = list(V = diag(1), nu = 0.002))

prior.flat.phylo <- list(G = list(
  G1 = list(V = diag(1)*0.01, nu = 0.01),
  G2 = list(V = diag(1)*0.01, nu = 0.01)),
  R = list(V = diag(1), nu = 0.002))

prior.iw.phylo <- list(G = list(
  G1 = list(V = diag(1)*0.002/2.002, nu = 2.002),
  G2 = list(V = diag(1)*0.002/2.002, nu = 2.002)),
  R = list(V = diag(1), nu = 0.002))



prior.exp.rand <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1),
            alpha.V = diag(1) * 1000),
  G2 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1),
            alpha.V = diag(1) * 1000),
  G3 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), 
            alpha.V = diag(1) * 1000)),
  R = list(V = diag(1), nu = 0.002))

prior.flat.rand <- list(G = list(
  G1 = list(V = diag(1)*0.01, nu = 0.01),
  G2 = list(V = diag(1)*0.01, nu = 0.01),
  G3 = list(V = diag(1)*0.01, nu = 0.01)),
  R = list(V = diag(1), nu = 0.002))

prior.iw.rand <- list(G = list(
  G1 = list(V = diag(1)*0.002/2.002, nu = 2.002),
  G2 = list(V = diag(1)*0.002/2.002, nu = 2.002),
  G3 = list(V = diag(1)*0.002/2.002, nu = 2.002)),
  R = list(V = diag(1), nu = 0.002))

#data$Social_category = factor(data$Social_category, levels=c('BifurDev', 'LinearDev', 'Solitary', 'Subsocial'))
#data$Social_category = factor(data$Social_category, levels=c('Solitary', 'Subsocial', 'LinearDev', 'BifurDev'))
#data$Social_category = factor(data$Social_category, levels=c('Subsocial', 'Solitary', 'LinearDev', 'BifurDev'))
#data$Social_category = factor(data$Social_category, levels=c('LinearDev', 'Solitary', 'Subsocial', 'BifurDev'))
#str(data$Social_category)
#order_levels = "SoSuLB"
#order_levels = "SuSoLB"
#order_levels = "LSoSuB"


#data$Food_category = factor(data$Food_category, levels=c('Litter', 'Wood-nesting', 'Wood', 'Fungus', 'Humus', 'Soil'))
#data$Food_category = factor(data$Food_category, levels=c('Wood-nesting', 'Litter', 'Wood', 'Fungus', 'Humus', 'Soil'))
#data$Food_category = factor(data$Food_category, levels=c('Wood', 'Litter', 'Wood-nesting', 'Fungus', 'Humus', 'Soil'))
#data$Food_category = factor(data$Food_category, levels=c('Fungus','Litter', 'Wood-nesting', 'Wood',  'Humus', 'Soil'))
#data$Food_category = factor(data$Food_category, levels=c('Humus','Litter', 'Wood-nesting', 'Wood', 'Fungus',  'Soil'))
data$Food_category = factor(data$Food_category, levels=c('Soil', 'Litter', 'Wood-nesting', 'Wood', 'Fungus', 'Humus'))

str(data$Food_category)
#order_levels = "LNWFHS"
#order_levels = "NLWFHS"
#order_levels = "WLNFHS"
#order_levels = "FLNWHS"
#order_levels = "HLNWFS"
order_levels = "SLNWFH"


#----------------------------------- Model -------------------------------------

k = as.numeric(args[1])
print(k)

if (k == 1){
  print("Exponential Model relax vs colony size with expanded priors and without random IDs.Rdata")
  model.kcsu.exp.phylo <- MCMCglmm(k ~ Colony_size_upper_range,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 13000000,
                                   thin = 5000, 
                                   burnin = 3000000, 
                                   pr = TRUE)

  print("saving...")
  save(model.kcsu.exp.phylo, file = paste0("Exponential Model relax vs colony size with expanded priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 2){
  print("Exponential Model relax vs colony size with flat priors and without random IDs.Rdata")
  model.kcsu.flat.phylo<- MCMCglmm(k ~ Colony_size_upper_range,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kcsu.flat.phylo, file = paste0("Exponential Model relax vs colony size with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 3){
  print("Exponential Model relax vs colony size with iw priors and without random IDs.Rdata")
  model.kcsu.iw.phylo<- MCMCglmm(k ~ Colony_size_upper_range,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kcsu.iw.phylo, file = paste0("Exponential Model relax vs colony size with iw priors and without random IDs.Rdata"))
  print("Done")
}



if (k == 4){
  print("Exponential Model relax vs colony size with exp priors and with random IDs.Rdata")
  model.kcsu.exp.rand<- MCMCglmm(k ~ Colony_size_upper_range,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kcsu.exp.rand, file = paste0("Exponential Model relax vs colony size with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 5){
  print("Exponential Model relax vs colony size with flat priors and with random IDs.Rdata")

  model.kcsu.flat.rand<- MCMCglmm(k ~ Colony_size_upper_range,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 13000000,
                            thin = 5000, 
                            burnin = 3000000, 
                            pr = TRUE)
  print("saving...")
  save(model.kcsu.flat.rand, file = paste0("Exponential Model relax vs colony size with flat priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 6){
  print("Exponential Model relax vs colony size with iw priors and with random IDs.Rdata")

  model.kcsu.iw.rand<- MCMCglmm(k ~ Colony_size_upper_range,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 13000000,
                          thin = 5000, 
                          burnin = 3000000, 
                          pr = TRUE)
  print("saving...")
  save(model.kcsu.iw.rand, file = paste0("Exponential Model relax vs colony size with iw priors and with random IDs.Rdata"))
  print("Done")
}





if (k == 7){
  print("Exponential Model relax vs ontogeny complexity metric with exp priors and without random IDs.Rdata")
  model.kocm.exp.phylo <- MCMCglmm(k ~ Ontogenic_complexity_metric,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 13000000,
                                   thin = 5000, 
                                   burnin = 3000000, 
                                   pr = TRUE)
  print("saving...")
  save(model.kocm.exp.phylo, file = paste0("Exponential Model relax vs ontogeny complexity metric with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 8){
  print("Exponential Model relax vs ontogeny complexity metric with flat priors and without random IDs.Rdata")

  model.kocm.flat.phylo<- MCMCglmm(k ~ Ontogenic_complexity_metric,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kocm.flat.phylo, file = paste0("Exponential Model relax vs ontogeny complexity metric with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 9){
  print("xponential Model relax vs ontogeny complexity metric with iw priors and without random IDs.Rdata")

  model.kocm.iw.phylo<- MCMCglmm(k ~ Ontogenic_complexity_metric,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kocm.iw.phylo, file = paste0("Exponential Model relax vs ontogeny complexity metric with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 10){
  print("Exponential Model relax vs ontogeny complexity metric with exp priors and with random IDs.Rdata")

  model.kocm.exp.rand<- MCMCglmm(k ~ Ontogenic_complexity_metric,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kocm.exp.rand, file = paste0("Exponential Model relax vs ontogeny complexity metric with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 11){
  print("Exponential Model relax vs ontogeny complexity metric with flat priors and with random IDs.Rdata")

  model.kocm.flat.rand<- MCMCglmm(k ~ Ontogenic_complexity_metric,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 13000000,
                            thin = 5000, 
                            burnin = 3000000, 
                            pr = TRUE)
  print("saving...")
  save(model.kocm.flat.rand, file = paste0("Exponential Model relax vs ontogeny complexity metric with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 12){
  print("Exponential Model relax vs ontogeny complexity metric with iw priors and with random IDs.Rdata")

  model.kocm.iw.rand<- MCMCglmm(k ~ Ontogenic_complexity_metric,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 13000000,
                          thin = 5000, 
                          burnin = 3000000, 
                          pr = TRUE)

  print("saving...")
  save(model.kocm.iw.rand, file = paste0("Exponential Model relax vs ontogeny complexity metric with iw priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 13){
  print("Exponential Model relax vs Ontogeny with expanded priors and without random IDs.Rdata")
  model.kont.exp.phylo <- MCMCglmm(k ~ Ontogeny,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 13000000,
                                   thin = 5000, 
                                   burnin = 3000000, 
                                   pr = TRUE)

  print("saving...")
  save(model.kont.exp.phylo, file = paste0("Exponential Model relax vs Ontogeny with expanded priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 14){
  print("Exponential Model relax vs Ontogeny with flat priors and without random IDs.Rdata")
  model.kont.flat.phylo<- MCMCglmm(k ~ Ontogeny,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kont.flat.phylo, file = paste0("Exponential Model relax vs Ontogeny with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 15){
  print("Exponential Model relax vs Ontogeny with iw priors and without random IDs.Rdata")
  model.kont.iw.phylo<- MCMCglmm(k ~ Ontogeny,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kont.iw.phylo, file = paste0("Exponential Model relax vs Ontogeny with iw priors and without random IDs.Rdata"))
  print("Done")
}



if (k == 16){
  print("Exponential Model relax vs Ontogeny with exp priors and with random IDs.Rdata")
  model.kont.exp.rand<- MCMCglmm(k ~ Ontogeny,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kont.exp.rand, file = paste0("Exponential Model relax vs Ontogeny with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 17){
  print("Exponential Model relax vs Ontogeny with flat priors and with random IDs.Rdata")

  model.kont.flat.rand<- MCMCglmm(k ~ Ontogeny,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 13000000,
                            thin = 5000, 
                            burnin = 3000000, 
                            pr = TRUE)
  print("saving...")
  save(model.kont.flat.rand, file = paste0("Exponential Model relax vs Ontogeny with flat priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 18){
  print("Exponential Model relax vs Ontogeny with iw priors and with random IDs.Rdata")

  model.kont.iw.rand<- MCMCglmm(k ~ Ontogeny,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 13000000,
                          thin = 5000, 
                          burnin = 3000000, 
                          pr = TRUE)
  print("saving...")
  save(model.kont.iw.rand, file = paste0("Exponential Model relax vs Ontogeny with iw priors and with random IDs.Rdata"))
  print("Done")
}





if (k == 19){
  print("Exponential Model relax vs number of terminal castes with exp priors and without random IDs.Rdata")
  model.kntc.exp.phylo <- MCMCglmm(k ~ No_terminal_caste,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 13000000,
                                   thin = 5000, 
                                   burnin = 3000000, 
                                   pr = TRUE)
  print("saving...")
  save(model.kntc.exp.phylo, file = paste0("Exponential Model relax vs number of terminal castes with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 20){
  print("Exponential Model relax vs number of terminal castes with flat priors and without random IDs.Rdata")

  model.kntc.flat.phylo<- MCMCglmm(k ~ No_terminal_caste,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kntc.flat.phylo, file = paste0("Exponential Model relax vs number of terminal castes with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 21){
  print("Exponential Model relax vs number of terminal castes with iw priors and without random IDs.Rdata")

  model.kntc.iw.phylo<- MCMCglmm(k ~ No_terminal_caste,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kntc.iw.phylo, file = paste0("Exponential Model relax vs number of terminal castes with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 22){
  print("Exponential Model relax vs number of terminal castes with exp priors and with random IDs.Rdata")

  model.kntc.exp.rand<- MCMCglmm(k ~ No_terminal_caste,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kntc.exp.rand, file = paste0("Exponential Model relax vs number of terminal castes with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 23){
  print("Exponential Model relax vs number of terminal castes with flat priors and with random IDs.Rdata")

  model.kntc.flat.rand<- MCMCglmm(k ~ No_terminal_caste,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 13000000,
                            thin = 5000, 
                            burnin = 3000000, 
                            pr = TRUE)
  print("saving...")
  save(model.kntc.flat.rand, file = paste0("Exponential Model relax vs number of terminal castes with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 24){
  print("Exponential Model relax vs number of terminal castes with iw priors and with random IDs.Rdata")

  model.kntc.iw.rand<- MCMCglmm(k ~ No_terminal_caste,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 13000000,
                          thin = 5000, 
                          burnin = 3000000, 
                          pr = TRUE)

  print("saving...")
  save(model.kntc.iw.rand, file = paste0("Exponential Model relax vs number of terminal castes with iw priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 25){
  print("Exponential Model relax vs social category with exp priors and without random IDs.Rdata")
  model.ksoc.exp.phylo <- MCMCglmm(k ~ Social_category,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 13000000,
                                   thin = 5000, 
                                   burnin = 3000000, 
                                   pr = TRUE)
  print("saving...")
  save(model.ksoc.exp.phylo, file = paste0("Exponential Model relax vs social category ", order_levels," with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 26){
  print("Exponential Model relax vs social category with flat priors and without random IDs.Rdata")

  model.ksoc.flat.phylo<- MCMCglmm(k ~ Social_category,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.ksoc.flat.phylo, file = paste0("Exponential Model relax vs social category ", order_levels," with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 27){
  print("Exponential Model relax vs social category with iw priors and without random IDs.Rdata")

  model.ksoc.iw.phylo<- MCMCglmm(k ~ Social_category,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.ksoc.iw.phylo, file = paste0("Exponential Model relax vs social category ", order_levels," with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 28){
  print("Exponential Model relax vs social category with exp priors and with random IDs.Rdata")

  model.ksoc.exp.rand<- MCMCglmm(k ~ Social_category,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.ksoc.exp.rand, file = paste0("Exponential Model relax vs social category ", order_levels," with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 29){
  print("Exponential Model relax vs social category with flat priors and with random IDs.Rdata")

  model.ksoc.flat.rand<- MCMCglmm(k ~ Social_category,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 13000000,
                            thin = 5000, 
                            burnin = 3000000, 
                            pr = TRUE)
  print("saving...")
  save(model.ksoc.flat.rand, file = paste0("Exponential Model relax vs social category ", order_levels," with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 30){
  print("Exponential Model relax vs social category with iw priors and with random IDs.Rdata")

  model.ksoc.iw.rand<- MCMCglmm(k ~ Social_category,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 13000000,
                          thin = 5000, 
                          burnin = 3000000, 
                          pr = TRUE)

  print("saving...")
  save(model.ksoc.iw.rand, file = paste0("Exponential Model relax vs social category ", order_levels," with iw priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 31){
  print("Exponential Model relax vs foraging with exp priors and without random IDs.Rdata")
  model.kfor.exp.phylo <- MCMCglmm(k ~ Foraging,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 13000000,
                                   thin = 5000, 
                                   burnin = 3000000, 
                                   pr = TRUE)
  print("saving...")
  save(model.kfor.exp.phylo, file = paste0("Exponential Model relax vs foraging with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 32){
  print("Exponential Model relax vs foraging with flat priors and without random IDs.Rdata")

  model.kfor.flat.phylo<- MCMCglmm(k ~ Foraging,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kfor.flat.phylo, file = paste0("Exponential Model relax vs foraging with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 33){
  print("Exponential Model relax vs foraging with iw priors and without random IDs.Rdata")

  model.kfor.iw.phylo<- MCMCglmm(k ~ Foraging,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kfor.iw.phylo, file = paste0("Exponential Model relax vs foraging with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 34){
  print("Exponential Model relax vs foraging with exp priors and with random IDs.Rdata")

  model.kfor.exp.rand<- MCMCglmm(k ~ Foraging,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kfor.exp.rand, file = paste0("Exponential Model relax vs foraging with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 35){
  print("Exponential Model relax vs foraging with flat priors and with random IDs.Rdata")

  model.kfor.flat.rand<- MCMCglmm(k ~ Foraging,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 13000000,
                            thin = 5000, 
                            burnin = 3000000, 
                            pr = TRUE)
  print("saving...")
  save(model.kfor.flat.rand, file = paste0("Exponential Model relax vs foraging with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 36){
  print("Exponential Model relax vs foraging with iw priors and with random IDs.Rdata")

  model.kfor.iw.rand<- MCMCglmm(k ~ Foraging,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 13000000,
                          thin = 5000, 
                          burnin = 3000000, 
                          pr = TRUE)

  print("saving...")
  save(model.kfor.iw.rand, file = paste0("Exponential Model relax vs foraging with iw priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 37){
  print("Exponential Model relax vs food category with exp priors and without random IDs.Rdata")
  model.kfca.exp.phylo <- MCMCglmm(k ~ Food_category,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 13000000,
                                   thin = 5000, 
                                   burnin = 3000000, 
                                   pr = TRUE)
  print("saving...")
  save(model.kfca.exp.phylo, file = paste0("Exponential Model relax vs food category ",order_levels," with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 38){
  print("Exponential Model relax vs food category with flat priors and without random IDs.Rdata")

  model.kfca.flat.phylo<- MCMCglmm(k ~ Food_category,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kfca.flat.phylo, file = paste0("Exponential Model relax vs food category ",order_levels," with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 39){
  print("Exponential Model relax vs food category with iw priors and without random IDs.Rdata")

  model.kfca.iw.phylo<- MCMCglmm(k ~ Food_category,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kfca.iw.phylo, file = paste0("Exponential Model relax vs food category ",order_levels," with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 40){
  print("Exponential Model relax vs food category with exp priors and with random IDs.Rdata")

  model.kfca.exp.rand<- MCMCglmm(k ~ Food_category,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kfca.exp.rand, file = paste0("Exponential Model relax vs food category ",order_levels," with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 41){
  print("Exponential Model relax vs food category with flat priors and with random IDs.Rdata")

  model.kfca.flat.rand<- MCMCglmm(k ~ Food_category,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 13000000,
                            thin = 5000, 
                            burnin = 3000000, 
                            pr = TRUE)
  print("saving...")
  save(model.kfca.flat.rand, file = paste0("Exponential Model relax vs food category ",order_levels," with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 42){
  print("Exponential Model relax vs food category with iw priors and with random IDs.Rdata")

  model.kfca.iw.rand<- MCMCglmm(k ~ Food_category,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 13000000,
                          thin = 5000, 
                          burnin = 3000000, 
                          pr = TRUE)

  print("saving...")
  save(model.kfca.iw.rand, file = paste0("Exponential Model relax vs food category ",order_levels," with iw priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 43){
  print("Exponential Model relax vs shortest path with exp priors and without random IDs.Rdata")
  model.kshp.exp.phylo <- MCMCglmm(k ~ Shortest_path,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 13000000,
                                   thin = 5000, 
                                   burnin = 3000000, 
                                   pr = TRUE)
  print("saving...")
  save(model.kshp.exp.phylo, file = paste0("Exponential Model relax vs shortest path with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 44){
  print("Exponential Model relax vs shortest path with flat priors and without random IDs.Rdata")

  model.kshp.flat.phylo<- MCMCglmm(k ~ Shortest_path,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kshp.flat.phylo, file = paste0("Exponential Model relax vs shortest path with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 45){
  print("Exponential Model relax vs shortest path with iw priors and without random IDs.Rdata")

  model.kshp.iw.phylo<- MCMCglmm(k ~ Shortest_path,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kshp.iw.phylo, file = paste0("Exponential Model relax vs shortest path with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 46){
  print("Exponential Model relax vs shortest path with exp priors and with random IDs.Rdata")

  model.kshp.exp.rand<- MCMCglmm(k ~ Shortest_path,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kshp.exp.rand, file = paste0("Exponential Model relax vs shortest path with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 47){
  print("Exponential Model relax vs shortest path with flat priors and with random IDs.Rdata")

  model.kshp.flat.rand<- MCMCglmm(k ~ Shortest_path,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 13000000,
                            thin = 5000, 
                            burnin = 3000000, 
                            pr = TRUE)
  print("saving...")
  save(model.kshp.flat.rand, file = paste0("Exponential Model relax vs shortest path with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 48){
  print("Exponential Model relax vs shortest path with iw priors and with random IDs.Rdata")

  model.kshp.iw.rand<- MCMCglmm(k ~ Shortest_path,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 13000000,
                          thin = 5000, 
                          burnin = 3000000, 
                          pr = TRUE)

  print("saving...")
  save(model.kshp.iw.rand, file = paste0("Exponential Model relax vs shortest path with iw priors and with random IDs.Rdata"))
  print("Done")
}



if (k == 49){
  print("Exponential Model relax vs protist with exp priors and without random IDs.Rdata")
  model.kpro.exp.phylo <- MCMCglmm(k ~ Protist,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 13000000,
                                   thin = 5000, 
                                   burnin = 3000000, 
                                   pr = TRUE)
  print("saving...")
  save(model.kpro.exp.phylo, file = paste0("Exponential Model relax vs protist with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 50){
  print("Exponential Model relax vs protist with flat priors and without random IDs.Rdata")

  model.kpro.flat.phylo<- MCMCglmm(k ~ Protist,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kpro.flat.phylo, file = paste0("Exponential Model relax vs protist with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 51){
  print("Exponential Model relax vs protist with iw priors and without random IDs.Rdata")

  model.kpro.iw.phylo<- MCMCglmm(k ~ Protist,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kpro.iw.phylo, file = paste0("Exponential Model relax vs protist with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 52){
  print("Exponential Model relax vs protist with exp priors and with random IDs.Rdata")

  model.kpro.exp.rand<- MCMCglmm(k ~ Protist,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 13000000,
                           thin = 5000, 
                           burnin = 3000000, 
                           pr = TRUE)
  print("saving...")
  save(model.kpro.exp.rand, file = paste0("Exponential Model relax vs protist with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 53){
  print("Exponential Model relax vs protist with flat priors and with random IDs.Rdata")

  model.kpro.flat.rand<- MCMCglmm(k ~ Protist,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 13000000,
                            thin = 5000, 
                            burnin = 3000000, 
                            pr = TRUE)
  print("saving...")
  save(model.kpro.flat.rand, file = paste0("Exponential Model relax vs protist with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 54){
  print("Exponential Model relax vs protist with iw priors and with random IDs.Rdata")

  model.kpro.iw.rand<- MCMCglmm(k ~ Protist,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 13000000,
                          thin = 5000, 
                          burnin = 3000000, 
                          pr = TRUE)

  print("saving...")
  save(model.kpro.iw.rand, file = paste0("Exponential Model relax vs protist with iw priors and with random IDs.Rdata"))
  print("Done")
}
