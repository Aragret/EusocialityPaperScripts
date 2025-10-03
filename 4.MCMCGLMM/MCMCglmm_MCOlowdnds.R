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

data = read.csv("data/dnds_MCO_mcmcglmm_dataset.csv",
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
  print("Exponential Model MCOdnds vs colony size with expanded priors and without random IDs.Rdata")
  model.mcsu.exp.phylo <- MCMCglmm(w_ratios_Value ~ Colony_size_upper_range,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 1300000,
                                   thin = 500, 
                                   burnin = 300000, 
                                   pr = TRUE)

  print("saving...")
  save(model.mcsu.exp.phylo, file = paste0("Exponential Model MCOdnds vs colony size with expanded priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 2){
  print("Exponential Model MCOdnds vs colony size with flat priors and without random IDs.Rdata")
  model.mcsu.flat.phylo<- MCMCglmm(w_ratios_Value ~ Colony_size_upper_range,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mcsu.flat.phylo, file = paste0("Exponential Model MCOdnds vs colony size with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 3){
  print("Exponential Model MCOdnds vs colony size with iw priors and without random IDs.Rdata")
  model.mcsu.iw.phylo<- MCMCglmm(w_ratios_Value ~ Colony_size_upper_range,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mcsu.iw.phylo, file = paste0("Exponential Model MCOdnds vs colony size with iw priors and without random IDs.Rdata"))
  print("Done")
}



if (k == 4){
  print("Exponential Model MCOdnds vs colony size with exp priors and with random IDs.Rdata")
  model.mcsu.exp.rand<- MCMCglmm(w_ratios_Value ~ Colony_size_upper_range,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mcsu.exp.rand, file = paste0("Exponential Model MCOdnds vs colony size with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 5){
  print("Exponential Model MCOdnds vs colony size with flat priors and with random IDs.Rdata")

  model.mcsu.flat.rand<- MCMCglmm(w_ratios_Value ~ Colony_size_upper_range,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 1300000,
                            thin = 500, 
                            burnin = 300000, 
                            pr = TRUE)
  print("saving...")
  save(model.mcsu.flat.rand, file = paste0("Exponential Model MCOdnds vs colony size with flat priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 6){
  print("Exponential Model MCOdnds vs colony size with iw priors and with random IDs.Rdata")

  model.mcsu.iw.rand<- MCMCglmm(w_ratios_Value ~ Colony_size_upper_range,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 1300000,
                          thin = 500, 
                          burnin = 300000, 
                          pr = TRUE)
  print("saving...")
  save(model.mcsu.iw.rand, file = paste0("Exponential Model MCOdnds vs colony size with iw priors and with random IDs.Rdata"))
  print("Done")
}





if (k == 7){
  print("Exponential Model MCOdnds vs ontogeny complexity metric with exp priors and without random IDs.Rdata")
  model.mocm.exp.phylo <- MCMCglmm(w_ratios_Value ~ Ontogenic_complexity_metric,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 1300000,
                                   thin = 500, 
                                   burnin = 300000, 
                                   pr = TRUE)
  print("saving...")
  save(model.mocm.exp.phylo, file = paste0("Exponential Model MCOdnds vs ontogeny complexity metric with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 8){
  print("Exponential Model MCOdnds vs ontogeny complexity metric with flat priors and without random IDs.Rdata")

  model.mocm.flat.phylo<- MCMCglmm(w_ratios_Value ~ Ontogenic_complexity_metric,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mocm.flat.phylo, file = paste0("Exponential Model MCOdnds vs ontogeny complexity metric with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 9){
  print("xponential Model MCOdnds vs ontogeny complexity metric with iw priors and without random IDs.Rdata")

  model.mocm.iw.phylo<- MCMCglmm(w_ratios_Value ~ Ontogenic_complexity_metric,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mocm.iw.phylo, file = paste0("Exponential Model MCOdnds vs ontogeny complexity metric with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 10){
  print("Exponential Model MCOdnds vs ontogeny complexity metric with exp priors and with random IDs.Rdata")

  model.mocm.exp.rand<- MCMCglmm(w_ratios_Value ~ Ontogenic_complexity_metric,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mocm.exp.rand, file = paste0("Exponential Model MCOdnds vs ontogeny complexity metric with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 11){
  print("Exponential Model MCOdnds vs ontogeny complexity metric with flat priors and with random IDs.Rdata")

  model.mocm.flat.rand<- MCMCglmm(w_ratios_Value ~ Ontogenic_complexity_metric,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 1300000,
                            thin = 500, 
                            burnin = 300000, 
                            pr = TRUE)
  print("saving...")
  save(model.mocm.flat.rand, file = paste0("Exponential Model MCOdnds vs ontogeny complexity metric with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 12){
  print("Exponential Model MCOdnds vs ontogeny complexity metric with iw priors and with random IDs.Rdata")

  model.mocm.iw.rand<- MCMCglmm(w_ratios_Value ~ Ontogenic_complexity_metric,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 1300000,
                          thin = 500, 
                          burnin = 300000, 
                          pr = TRUE)

  print("saving...")
  save(model.mocm.iw.rand, file = paste0("Exponential Model MCOdnds vs ontogeny complexity metric with iw priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 13){
  print("Exponential Model MCOdnds vs Ontogeny with expanded priors and without random IDs.Rdata")
  model.mont.exp.phylo <- MCMCglmm(w_ratios_Value ~ Ontogeny,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 1300000,
                                   thin = 500, 
                                   burnin = 300000, 
                                   pr = TRUE)

  print("saving...")
  save(model.mont.exp.phylo, file = paste0("Exponential Model MCOdnds vs Ontogeny with expanded priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 14){
  print("Exponential Model MCOdnds vs Ontogeny with flat priors and without random IDs.Rdata")
  model.mont.flat.phylo<- MCMCglmm(w_ratios_Value ~ Ontogeny,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mont.flat.phylo, file = paste0("Exponential Model MCOdnds vs Ontogeny with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 15){
  print("Exponential Model MCOdnds vs Ontogeny with iw priors and without random IDs.Rdata")
  model.mont.iw.phylo<- MCMCglmm(w_ratios_Value ~ Ontogeny,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mont.iw.phylo, file = paste0("Exponential Model MCOdnds vs Ontogeny with iw priors and without random IDs.Rdata"))
  print("Done")
}



if (k == 16){
  print("Exponential Model MCOdnds vs Ontogeny with exp priors and with random IDs.Rdata")
  model.mont.exp.rand<- MCMCglmm(w_ratios_Value ~ Ontogeny,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mont.exp.rand, file = paste0("Exponential Model MCOdnds vs Ontogeny with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 17){
  print("Exponential Model MCOdnds vs Ontogeny with flat priors and with random IDs.Rdata")

  model.mont.flat.rand<- MCMCglmm(w_ratios_Value ~ Ontogeny,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 1300000,
                            thin = 500, 
                            burnin = 300000, 
                            pr = TRUE)
  print("saving...")
  save(model.mont.flat.rand, file = paste0("Exponential Model MCOdnds vs Ontogeny with flat priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 18){
  print("Exponential Model MCOdnds vs Ontogeny with iw priors and with random IDs.Rdata")

  model.mont.iw.rand<- MCMCglmm(w_ratios_Value ~ Ontogeny,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 1300000,
                          thin = 500, 
                          burnin = 300000, 
                          pr = TRUE)
  print("saving...")
  save(model.mont.iw.rand, file = paste0("Exponential Model MCOdnds vs Ontogeny with iw priors and with random IDs.Rdata"))
  print("Done")
}





if (k == 19){
  print("Exponential Model MCOdnds vs number of terminal castes with exp priors and without random IDs.Rdata")
  model.mntc.exp.phylo <- MCMCglmm(w_ratios_Value ~ No_terminal_caste,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 1300000,
                                   thin = 500, 
                                   burnin = 300000, 
                                   pr = TRUE)
  print("saving...")
  save(model.mntc.exp.phylo, file = paste0("Exponential Model MCOdnds vs number of terminal castes with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 20){
  print("Exponential Model MCOdnds vs number of terminal castes with flat priors and without random IDs.Rdata")

  model.mntc.flat.phylo<- MCMCglmm(w_ratios_Value ~ No_terminal_caste,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mntc.flat.phylo, file = paste0("Exponential Model MCOdnds vs number of terminal castes with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 21){
  print("Exponential Model MCOdnds vs number of terminal castes with iw priors and without random IDs.Rdata")

  model.mntc.iw.phylo<- MCMCglmm(w_ratios_Value ~ No_terminal_caste,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mntc.iw.phylo, file = paste0("Exponential Model MCOdnds vs number of terminal castes with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 22){
  print("Exponential Model MCOdnds vs number of terminal castes with exp priors and with random IDs.Rdata")

  model.mntc.exp.rand<- MCMCglmm(w_ratios_Value ~ No_terminal_caste,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mntc.exp.rand, file = paste0("Exponential Model MCOdnds vs number of terminal castes with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 23){
  print("Exponential Model MCOdnds vs number of terminal castes with flat priors and with random IDs.Rdata")

  model.mntc.flat.rand<- MCMCglmm(w_ratios_Value ~ No_terminal_caste,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 1300000,
                            thin = 500, 
                            burnin = 300000, 
                            pr = TRUE)
  print("saving...")
  save(model.mntc.flat.rand, file = paste0("Exponential Model MCOdnds vs number of terminal castes with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 24){
  print("Exponential Model MCOdnds vs number of terminal castes with iw priors and with random IDs.Rdata")

  model.mntc.iw.rand<- MCMCglmm(w_ratios_Value ~ No_terminal_caste,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 1300000,
                          thin = 500, 
                          burnin = 300000, 
                          pr = TRUE)

  print("saving...")
  save(model.mntc.iw.rand, file = paste0("Exponential Model MCOdnds vs number of terminal castes with iw priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 25){
  print("Exponential Model MCOdnds vs social category with exp priors and without random IDs.Rdata")
  model.msoc.exp.phylo <- MCMCglmm(w_ratios_Value ~ Social_category,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 1300000,
                                   thin = 500, 
                                   burnin = 300000, 
                                   pr = TRUE)
  print("saving...")
  save(model.msoc.exp.phylo, file = paste0("Exponential Model MCOdnds vs social category ", order_levels," with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 26){
  print("Exponential Model MCOdnds vs social category with flat priors and without random IDs.Rdata")

  model.msoc.flat.phylo<- MCMCglmm(w_ratios_Value ~ Social_category,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.msoc.flat.phylo, file = paste0("Exponential Model MCOdnds vs social category ", order_levels," with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 27){
  print("Exponential Model MCOdnds vs social category with iw priors and without random IDs.Rdata")

  model.msoc.iw.phylo<- MCMCglmm(w_ratios_Value ~ Social_category,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.msoc.iw.phylo, file = paste0("Exponential Model MCOdnds vs social category ", order_levels," with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 28){
  print("Exponential Model MCOdnds vs social category with exp priors and with random IDs.Rdata")

  model.msoc.exp.rand<- MCMCglmm(w_ratios_Value ~ Social_category,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.msoc.exp.rand, file = paste0("Exponential Model MCOdnds vs social category ", order_levels," with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 29){
  print("Exponential Model MCOdnds vs social category with flat priors and with random IDs.Rdata")

  model.msoc.flat.rand<- MCMCglmm(w_ratios_Value ~ Social_category,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 1300000,
                            thin = 500, 
                            burnin = 300000, 
                            pr = TRUE)
  print("saving...")
  save(model.msoc.flat.rand, file = paste0("Exponential Model MCOdnds vs social category ", order_levels," with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 30){
  print("Exponential Model MCOdnds vs social category with iw priors and with random IDs.Rdata")

  model.msoc.iw.rand<- MCMCglmm(w_ratios_Value ~ Social_category,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 1300000,
                          thin = 500, 
                          burnin = 300000, 
                          pr = TRUE)

  print("saving...")
  save(model.msoc.iw.rand, file = paste0("Exponential Model MCOdnds vs social category ", order_levels," with iw priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 31){
  print("Exponential Model MCOdnds vs foraging with exp priors and without random IDs.Rdata")
  model.mfor.exp.phylo <- MCMCglmm(w_ratios_Value ~ Foraging,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 1300000,
                                   thin = 500, 
                                   burnin = 300000, 
                                   pr = TRUE)
  print("saving...")
  save(model.mfor.exp.phylo, file = paste0("Exponential Model MCOdnds vs foraging with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 32){
  print("Exponential Model MCOdnds vs foraging with flat priors and without random IDs.Rdata")

  model.mfor.flat.phylo<- MCMCglmm(w_ratios_Value ~ Foraging,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mfor.flat.phylo, file = paste0("Exponential Model MCOdnds vs foraging with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 33){
  print("Exponential Model MCOdnds vs foraging with iw priors and without random IDs.Rdata")

  model.mfor.iw.phylo<- MCMCglmm(w_ratios_Value ~ Foraging,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mfor.iw.phylo, file = paste0("Exponential Model MCOdnds vs foraging with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 34){
  print("Exponential Model MCOdnds vs foraging with exp priors and with random IDs.Rdata")

  model.mfor.exp.rand<- MCMCglmm(w_ratios_Value ~ Foraging,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mfor.exp.rand, file = paste0("Exponential Model MCOdnds vs foraging with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 35){
  print("Exponential Model MCOdnds vs foraging with flat priors and with random IDs.Rdata")

  model.mfor.flat.rand<- MCMCglmm(w_ratios_Value ~ Foraging,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 1300000,
                            thin = 500, 
                            burnin = 300000, 
                            pr = TRUE)
  print("saving...")
  save(model.mfor.flat.rand, file = paste0("Exponential Model MCOdnds vs foraging with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 36){
  print("Exponential Model MCOdnds vs foraging with iw priors and with random IDs.Rdata")

  model.mfor.iw.rand<- MCMCglmm(w_ratios_Value ~ Foraging,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 1300000,
                          thin = 500, 
                          burnin = 300000, 
                          pr = TRUE)

  print("saving...")
  save(model.mfor.iw.rand, file = paste0("Exponential Model MCOdnds vs foraging with iw priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 37){
  print("Exponential Model MCOdnds vs food category with exp priors and without random IDs.Rdata")
  model.mfca.exp.phylo <- MCMCglmm(w_ratios_Value ~ Food_category,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 1300000,
                                   thin = 500, 
                                   burnin = 300000, 
                                   pr = TRUE)
  print("saving...")
  save(model.mfca.exp.phylo, file = paste0("Exponential Model MCOdnds vs food category ",order_levels," with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 38){
  print("Exponential Model MCOdnds vs food category with flat priors and without random IDs.Rdata")

  model.mfca.flat.phylo<- MCMCglmm(w_ratios_Value ~ Food_category,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mfca.flat.phylo, file = paste0("Exponential Model MCOdnds vs food category ",order_levels," with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 39){
  print("Exponential Model MCOdnds vs food category with iw priors and without random IDs.Rdata")

  model.mfca.iw.phylo<- MCMCglmm(w_ratios_Value ~ Food_category,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mfca.iw.phylo, file = paste0("Exponential Model MCOdnds vs food category ",order_levels," with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 40){
  print("Exponential Model MCOdnds vs food category with exp priors and with random IDs.Rdata")

  model.mfca.exp.rand<- MCMCglmm(w_ratios_Value ~ Food_category,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mfca.exp.rand, file = paste0("Exponential Model MCOdnds vs food category ",order_levels," with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 41){
  print("Exponential Model MCOdnds vs food category with flat priors and with random IDs.Rdata")

  model.mfca.flat.rand<- MCMCglmm(w_ratios_Value ~ Food_category,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 1300000,
                            thin = 500, 
                            burnin = 300000, 
                            pr = TRUE)
  print("saving...")
  save(model.mfca.flat.rand, file = paste0("Exponential Model MCOdnds vs food category ",order_levels," with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 42){
  print("Exponential Model MCOdnds vs food category with iw priors and with random IDs.Rdata")

  model.mfca.iw.rand<- MCMCglmm(w_ratios_Value ~ Food_category,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 1300000,
                          thin = 500, 
                          burnin = 300000, 
                          pr = TRUE)

  print("saving...")
  save(model.mfca.iw.rand, file = paste0("Exponential Model MCOdnds vs food category ",order_levels," with iw priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 43){
  print("Exponential Model MCOdnds vs shortest path with exp priors and without random IDs.Rdata")
  model.mshp.exp.phylo <- MCMCglmm(w_ratios_Value ~ Shortest_path,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 1300000,
                                   thin = 500, 
                                   burnin = 300000, 
                                   pr = TRUE)
  print("saving...")
  save(model.mshp.exp.phylo, file = paste0("Exponential Model MCOdnds vs shortest path with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 44){
  print("Exponential Model MCOdnds vs shortest path with flat priors and without random IDs.Rdata")

  model.mshp.flat.phylo<- MCMCglmm(w_ratios_Value ~ Shortest_path,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mshp.flat.phylo, file = paste0("Exponential Model MCOdnds vs shortest path with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 45){
  print("Exponential Model MCOdnds vs shortest path with iw priors and without random IDs.Rdata")

  model.mshp.iw.phylo<- MCMCglmm(w_ratios_Value ~ Shortest_path,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mshp.iw.phylo, file = paste0("Exponential Model MCOdnds vs shortest path with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 46){
  print("Exponential Model MCOdnds vs shortest path with exp priors and with random IDs.Rdata")

  model.mshp.exp.rand<- MCMCglmm(w_ratios_Value ~ Shortest_path,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mshp.exp.rand, file = paste0("Exponential Model MCOdnds vs shortest path with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 47){
  print("Exponential Model MCOdnds vs shortest path with flat priors and with random IDs.Rdata")

  model.mshp.flat.rand<- MCMCglmm(w_ratios_Value ~ Shortest_path,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 1300000,
                            thin = 500, 
                            burnin = 300000, 
                            pr = TRUE)
  print("saving...")
  save(model.mshp.flat.rand, file = paste0("Exponential Model MCOdnds vs shortest path with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 48){
  print("Exponential Model MCOdnds vs shortest path with iw priors and with random IDs.Rdata")

  model.mshp.iw.rand<- MCMCglmm(w_ratios_Value ~ Shortest_path,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 1300000,
                          thin = 500, 
                          burnin = 300000, 
                          pr = TRUE)

  print("saving...")
  save(model.mshp.iw.rand, file = paste0("Exponential Model MCOdnds vs shortest path with iw priors and with random IDs.Rdata"))
  print("Done")
}



if (k == 49){
  print("Exponential Model MCOdnds vs protist with exp priors and without random IDs.Rdata")
  model.mpro.exp.phylo <- MCMCglmm(w_ratios_Value ~ Protist,
                                   random = ~Species + Orthogroup,
                                   family = "exponential",
                                   ginverse = list(Species = inv.phylo), 
                                   prior = prior.exp.phylo,
                                   data = data,
                                   nitt = 1300000,
                                   thin = 500, 
                                   burnin = 300000, 
                                   pr = TRUE)
  print("saving...")
  save(model.mpro.exp.phylo, file = paste0("Exponential Model MCOdnds vs protist with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 50){
  print("Exponential Model MCOdnds vs protist with flat priors and without random IDs.Rdata")

  model.mpro.flat.phylo<- MCMCglmm(w_ratios_Value ~ Protist,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.flat.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mpro.flat.phylo, file = paste0("Exponential Model MCOdnds vs protist with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 51){
  print("Exponential Model MCOdnds vs protist with iw priors and without random IDs.Rdata")

  model.mpro.iw.phylo<- MCMCglmm(w_ratios_Value ~ Protist,
                           random = ~Species + Orthogroup,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.iw.phylo,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mpro.iw.phylo, file = paste0("Exponential Model MCOdnds vs protist with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 52){
  print("Exponential Model MCOdnds vs protist with exp priors and with random IDs.Rdata")

  model.mpro.exp.rand<- MCMCglmm(w_ratios_Value ~ Protist,
                           random = ~Species + Orthogroup + random_id,
                           family = "exponential",
                           ginverse = list(Species = inv.phylo), 
                           prior = prior.exp.rand,
                           data = data,
                           nitt = 1300000,
                           thin = 500, 
                           burnin = 300000, 
                           pr = TRUE)
  print("saving...")
  save(model.mpro.exp.rand, file = paste0("Exponential Model MCOdnds vs protist with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 53){
  print("Exponential Model MCOdnds vs protist with flat priors and with random IDs.Rdata")

  model.mpro.flat.rand<- MCMCglmm(w_ratios_Value ~ Protist,
                            random = ~Species + Orthogroup + random_id,
                            family = "exponential",
                            ginverse = list(Species = inv.phylo), 
                            prior = prior.flat.rand,
                            data = data,
                            nitt = 1300000,
                            thin = 500, 
                            burnin = 300000, 
                            pr = TRUE)
  print("saving...")
  save(model.mpro.flat.rand, file = paste0("Exponential Model MCOdnds vs protist with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 54){
  print("Exponential Model MCOdnds vs protist with iw priors and with random IDs.Rdata")

  model.mpro.iw.rand<- MCMCglmm(w_ratios_Value ~ Protist,
                          random = ~Species + Orthogroup + random_id,
                          family = "exponential",
                          ginverse = list(Species = inv.phylo), 
                          prior = prior.iw.rand,
                          data = data,
                          nitt = 1300000,
                          thin = 500, 
                          burnin = 300000, 
                          pr = TRUE)

  print("saving...")
  save(model.mpro.iw.rand, file = paste0("Exponential Model MCOdnds vs protist with iw priors and with random IDs.Rdata"))
  print("Done")
}

