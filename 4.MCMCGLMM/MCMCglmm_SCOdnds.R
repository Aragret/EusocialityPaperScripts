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

data = read.csv("data/dnds_SCO_mcmcglmm_dataset.csv",
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
  print("Exponential Model dnds vs colony size with expanded priors and without random IDs.Rdata")
  model.wcsu.exp.phylo <- MCMCglmm(w_ratios_Value ~ Colony_size_upper_range,
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
  save(model.wcsu.exp.phylo, file = paste0("Exponential Model dnds vs colony size with expanded priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 2){
  print("Exponential Model dnds vs colony size with flat priors and without random IDs.Rdata")
  model.wcsu.flat.phylo<- MCMCglmm(w_ratios_Value ~ Colony_size_upper_range,
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
  save(model.wcsu.flat.phylo, file = paste0("Exponential Model dnds vs colony size with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 3){
  print("Exponential Model dnds vs colony size with iw priors and without random IDs.Rdata")
  model.wcsu.iw.phylo<- MCMCglmm(w_ratios_Value ~ Colony_size_upper_range,
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
  save(model.wcsu.iw.phylo, file = paste0("Exponential Model dnds vs colony size with iw priors and without random IDs.Rdata"))
  print("Done")
}



if (k == 4){
  print("Exponential Model dnds vs colony size with exp priors and with random IDs.Rdata")
  model.wcsu.exp.rand<- MCMCglmm(w_ratios_Value ~ Colony_size_upper_range,
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
  save(model.wcsu.exp.rand, file = paste0("Exponential Model dnds vs colony size with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 5){
  print("Exponential Model dnds vs colony size with flat priors and with random IDs.Rdata")

  model.wcsu.flat.rand<- MCMCglmm(w_ratios_Value ~ Colony_size_upper_range,
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
  save(model.wcsu.flat.rand, file = paste0("Exponential Model dnds vs colony size with flat priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 6){
  print("Exponential Model dnds vs colony size with iw priors and with random IDs.Rdata")

  model.wcsu.iw.rand<- MCMCglmm(w_ratios_Value ~ Colony_size_upper_range,
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
  save(model.wcsu.iw.rand, file = paste0("Exponential Model dnds vs colony size with iw priors and with random IDs.Rdata"))
  print("Done")
}





if (k == 7){
  print("Exponential Model dnds vs ontogeny complexity metric with exp priors and without random IDs.Rdata")
  model.wocm.exp.phylo <- MCMCglmm(w_ratios_Value ~ Ontogenic_complexity_metric,
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
  save(model.wocm.exp.phylo, file = paste0("Exponential Model dnds vs ontogeny complexity metric with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 8){
  print("Exponential Model dnds vs ontogeny complexity metric with flat priors and without random IDs.Rdata")

  model.wocm.flat.phylo<- MCMCglmm(w_ratios_Value ~ Ontogenic_complexity_metric,
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
  save(model.wocm.flat.phylo, file = paste0("Exponential Model dnds vs ontogeny complexity metric with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 9){
  print("xponential Model dnds vs ontogeny complexity metric with iw priors and without random IDs.Rdata")

  model.wocm.iw.phylo<- MCMCglmm(w_ratios_Value ~ Ontogenic_complexity_metric,
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
  save(model.wocm.iw.phylo, file = paste0("Exponential Model dnds vs ontogeny complexity metric with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 10){
  print("Exponential Model dnds vs ontogeny complexity metric with exp priors and with random IDs.Rdata")

  model.wocm.exp.rand<- MCMCglmm(w_ratios_Value ~ Ontogenic_complexity_metric,
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
  save(model.wocm.exp.rand, file = paste0("Exponential Model dnds vs ontogeny complexity metric with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 11){
  print("Exponential Model dnds vs ontogeny complexity metric with flat priors and with random IDs.Rdata")

  model.wocm.flat.rand<- MCMCglmm(w_ratios_Value ~ Ontogenic_complexity_metric,
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
  save(model.wocm.flat.rand, file = paste0("Exponential Model dnds vs ontogeny complexity metric with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 12){
  print("Exponential Model dnds vs ontogeny complexity metric with iw priors and with random IDs.Rdata")

  model.wocm.iw.rand<- MCMCglmm(w_ratios_Value ~ Ontogenic_complexity_metric,
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
  save(model.wocm.iw.rand, file = paste0("Exponential Model dnds vs ontogeny complexity metric with iw priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 13){
  print("Exponential Model dnds vs Ontogeny with expanded priors and without random IDs.Rdata")
  model.wont.exp.phylo <- MCMCglmm(w_ratios_Value ~ Ontogeny,
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
  save(model.wont.exp.phylo, file = paste0("Exponential Model dnds vs Ontogeny with expanded priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 14){
  print("Exponential Model dnds vs Ontogeny with flat priors and without random IDs.Rdata")
  model.wont.flat.phylo<- MCMCglmm(w_ratios_Value ~ Ontogeny,
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
  save(model.wont.flat.phylo, file = paste0("Exponential Model dnds vs Ontogeny with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 15){
  print("Exponential Model dnds vs Ontogeny with iw priors and without random IDs.Rdata")
  model.wont.iw.phylo<- MCMCglmm(w_ratios_Value ~ Ontogeny,
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
  save(model.wont.iw.phylo, file = paste0("Exponential Model dnds vs Ontogeny with iw priors and without random IDs.Rdata"))
  print("Done")
}



if (k == 16){
  print("Exponential Model dnds vs Ontogeny with exp priors and with random IDs.Rdata")
  model.wont.exp.rand<- MCMCglmm(w_ratios_Value ~ Ontogeny,
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
  save(model.wont.exp.rand, file = paste0("Exponential Model dnds vs Ontogeny with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 17){
  print("Exponential Model dnds vs Ontogeny with flat priors and with random IDs.Rdata")

  model.wont.flat.rand<- MCMCglmm(w_ratios_Value ~ Ontogeny,
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
  save(model.wont.flat.rand, file = paste0("Exponential Model dnds vs Ontogeny with flat priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 18){
  print("Exponential Model dnds vs Ontogeny with iw priors and with random IDs.Rdata")

  model.wont.iw.rand<- MCMCglmm(w_ratios_Value ~ Ontogeny,
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
  save(model.wont.iw.rand, file = paste0("Exponential Model dnds vs Ontogeny with iw priors and with random IDs.Rdata"))
  print("Done")
}





if (k == 19){
  print("Exponential Model dnds vs number of terminal castes with exp priors and without random IDs.Rdata")
  model.wntc.exp.phylo <- MCMCglmm(w_ratios_Value ~ No_terminal_caste,
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
  save(model.wntc.exp.phylo, file = paste0("Exponential Model dnds vs number of terminal castes with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 20){
  print("Exponential Model dnds vs number of terminal castes with flat priors and without random IDs.Rdata")

  model.wntc.flat.phylo<- MCMCglmm(w_ratios_Value ~ No_terminal_caste,
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
  save(model.wntc.flat.phylo, file = paste0("Exponential Model dnds vs number of terminal castes with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 21){
  print("Exponential Model dnds vs number of terminal castes with iw priors and without random IDs.Rdata")

  model.wntc.iw.phylo<- MCMCglmm(w_ratios_Value ~ No_terminal_caste,
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
  save(model.wntc.iw.phylo, file = paste0("Exponential Model dnds vs number of terminal castes with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 22){
  print("Exponential Model dnds vs number of terminal castes with exp priors and with random IDs.Rdata")

  model.wntc.exp.rand<- MCMCglmm(w_ratios_Value ~ No_terminal_caste,
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
  save(model.wntc.exp.rand, file = paste0("Exponential Model dnds vs number of terminal castes with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 23){
  print("Exponential Model dnds vs number of terminal castes with flat priors and with random IDs.Rdata")

  model.wntc.flat.rand<- MCMCglmm(w_ratios_Value ~ No_terminal_caste,
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
  save(model.wntc.flat.rand, file = paste0("Exponential Model dnds vs number of terminal castes with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 24){
  print("Exponential Model dnds vs number of terminal castes with iw priors and with random IDs.Rdata")

  model.wntc.iw.rand<- MCMCglmm(w_ratios_Value ~ No_terminal_caste,
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
  save(model.wntc.iw.rand, file = paste0("Exponential Model dnds vs number of terminal castes with iw priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 25){
  print("Exponential Model dnds vs social category with exp priors and without random IDs.Rdata")
  model.wsoc.exp.phylo <- MCMCglmm(w_ratios_Value ~ Social_category,
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
  save(model.wsoc.exp.phylo, file = paste0("Exponential Model dnds vs social category ", order_levels," with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 26){
  print("Exponential Model dnds vs social category with flat priors and without random IDs.Rdata")

  model.wsoc.flat.phylo<- MCMCglmm(w_ratios_Value ~ Social_category,
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
  save(model.wsoc.flat.phylo, file = paste0("Exponential Model dnds vs social category ", order_levels," with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 27){
  print("Exponential Model dnds vs social category with iw priors and without random IDs.Rdata")

  model.wsoc.iw.phylo<- MCMCglmm(w_ratios_Value ~ Social_category,
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
  save(model.wsoc.iw.phylo, file = paste0("Exponential Model dnds vs social category ", order_levels," with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 28){
  print("Exponential Model dnds vs social category with exp priors and with random IDs.Rdata")

  model.wsoc.exp.rand<- MCMCglmm(w_ratios_Value ~ Social_category,
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
  save(model.wsoc.exp.rand, file = paste0("Exponential Model dnds vs social category ", order_levels," with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 29){
  print("Exponential Model dnds vs social category with flat priors and with random IDs.Rdata")

  model.wsoc.flat.rand<- MCMCglmm(w_ratios_Value ~ Social_category,
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
  save(model.wsoc.flat.rand, file = paste0("Exponential Model dnds vs social category ", order_levels," with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 30){
  print("Exponential Model dnds vs social category with iw priors and with random IDs.Rdata")

  model.wsoc.iw.rand<- MCMCglmm(w_ratios_Value ~ Social_category,
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
  save(model.wsoc.iw.rand, file = paste0("Exponential Model dnds vs social category ", order_levels," with iw priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 31){
  print("Exponential Model dnds vs foraging with exp priors and without random IDs.Rdata")
  model.wfor.exp.phylo <- MCMCglmm(w_ratios_Value ~ Foraging,
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
  save(model.wfor.exp.phylo, file = paste0("Exponential Model dnds vs foraging with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 32){
  print("Exponential Model dnds vs foraging with flat priors and without random IDs.Rdata")

  model.wfor.flat.phylo<- MCMCglmm(w_ratios_Value ~ Foraging,
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
  save(model.wfor.flat.phylo, file = paste0("Exponential Model dnds vs foraging with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 33){
  print("Exponential Model dnds vs foraging with iw priors and without random IDs.Rdata")

  model.wfor.iw.phylo<- MCMCglmm(w_ratios_Value ~ Foraging,
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
  save(model.wfor.iw.phylo, file = paste0("Exponential Model dnds vs foraging with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 34){
  print("Exponential Model dnds vs foraging with exp priors and with random IDs.Rdata")

  model.wfor.exp.rand<- MCMCglmm(w_ratios_Value ~ Foraging,
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
  save(model.wfor.exp.rand, file = paste0("Exponential Model dnds vs foraging with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 35){
  print("Exponential Model dnds vs foraging with flat priors and with random IDs.Rdata")

  model.wfor.flat.rand<- MCMCglmm(w_ratios_Value ~ Foraging,
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
  save(model.wfor.flat.rand, file = paste0("Exponential Model dnds vs foraging with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 36){
  print("Exponential Model dnds vs foraging with iw priors and with random IDs.Rdata")

  model.wfor.iw.rand<- MCMCglmm(w_ratios_Value ~ Foraging,
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
  save(model.wfor.iw.rand, file = paste0("Exponential Model dnds vs foraging with iw priors and with random IDs.Rdata"))
  print("Done")
}


if (k == 37){
  print("Exponential Model dnds vs food category with exp priors and without random IDs.Rdata")
  model.wfca.exp.phylo <- MCMCglmm(w_ratios_Value ~ Food_category,
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
  save(model.wfca.exp.phylo, file = paste0("Exponential Model dnds vs food category ", order_levels," with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 38){
  print("Exponential Model dnds vs food category with flat priors and without random IDs.Rdata")

  model.wfca.flat.phylo<- MCMCglmm(w_ratios_Value ~ Food_category,
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
  save(model.wfca.flat.phylo, file = paste0("Exponential Model dnds vs food category ", order_levels," with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 39){
  print("Exponential Model dnds vs food category with iw priors and without random IDs.Rdata")

  model.wfca.iw.phylo<- MCMCglmm(w_ratios_Value ~ Food_category,
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
  save(model.wfca.iw.phylo, file = paste0("Exponential Model dnds vs food category ", order_levels," with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 40){
  print("Exponential Model dnds vs food category with exp priors and with random IDs.Rdata")

  model.wfca.exp.rand<- MCMCglmm(w_ratios_Value ~ Food_category,
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
  save(model.wfca.exp.rand, file = paste0("Exponential Model dnds vs food category ", order_levels," with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 41){
  print("Exponential Model dnds vs food category with flat priors and with random IDs.Rdata")

  model.wfca.flat.rand<- MCMCglmm(w_ratios_Value ~ Food_category,
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
  save(model.wfca.flat.rand, file = paste0("Exponential Model dnds vs food category ", order_levels," with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 42){
  print("Exponential Model dnds vs food category with iw priors and with random IDs.Rdata")

  model.wfca.iw.rand<- MCMCglmm(w_ratios_Value ~ Food_category,
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
  save(model.wfca.iw.rand, file = paste0("Exponential Model dnds vs food category ", order_levels," with iw priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 43){
  print("Exponential Model dnds vs shortest path with exp priors and without random IDs.Rdata")
  model.wshp.exp.phylo <- MCMCglmm(w_ratios_Value ~ Shortest_path,
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
  save(model.wshp.exp.phylo, file = paste0("Exponential Model dnds vs shortest path with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 44){
  print("Exponential Model dnds vs shortest path with flat priors and without random IDs.Rdata")

  model.wshp.flat.phylo<- MCMCglmm(w_ratios_Value ~ Shortest_path,
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
  save(model.wshp.flat.phylo, file = paste0("Exponential Model dnds vs shortest path with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 45){
  print("Exponential Model dnds vs shortest path with iw priors and without random IDs.Rdata")

  model.wshp.iw.phylo<- MCMCglmm(w_ratios_Value ~ Shortest_path,
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
  save(model.wshp.iw.phylo, file = paste0("Exponential Model dnds vs shortest path with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 46){
  print("Exponential Model dnds vs shortest path with exp priors and with random IDs.Rdata")

  model.wshp.exp.rand<- MCMCglmm(w_ratios_Value ~ Shortest_path,
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
  save(model.wshp.exp.rand, file = paste0("Exponential Model dnds vs shortest path with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 47){
  print("Exponential Model dnds vs shortest path with flat priors and with random IDs.Rdata")

  model.wshp.flat.rand<- MCMCglmm(w_ratios_Value ~ Shortest_path,
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
  save(model.wshp.flat.rand, file = paste0("Exponential Model dnds vs shortest path with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 48){
  print("Exponential Model dnds vs shortest path with iw priors and with random IDs.Rdata")

  model.wshp.iw.rand<- MCMCglmm(w_ratios_Value ~ Shortest_path,
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
  save(model.wshp.iw.rand, file = paste0("Exponential Model dnds vs shortest path with iw priors and with random IDs.Rdata"))
  print("Done")
}



if (k == 49){
  print("Exponential Model dnds vs protist with exp priors and without random IDs.Rdata")
  model.wpro.exp.phylo <- MCMCglmm(w_ratios_Value ~ Protist,
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
  save(model.wpro.exp.phylo, file = paste0("Exponential Model dnds vs protist with exp priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 50){
  print("Exponential Model dnds vs protist with flat priors and without random IDs.Rdata")

  model.wpro.flat.phylo<- MCMCglmm(w_ratios_Value ~ Protist,
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
  save(model.wpro.flat.phylo, file = paste0("Exponential Model dnds vs protist with flat priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 51){
  print("Exponential Model dnds vs protist with iw priors and without random IDs.Rdata")

  model.wpro.iw.phylo<- MCMCglmm(w_ratios_Value ~ Protist,
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
  save(model.wpro.iw.phylo, file = paste0("Exponential Model dnds vs protist with iw priors and without random IDs.Rdata"))
  print("Done")
}


if (k == 52){
  print("Exponential Model dnds vs protist with exp priors and with random IDs.Rdata")

  model.wpro.exp.rand<- MCMCglmm(w_ratios_Value ~ Protist,
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
  save(model.wpro.exp.rand, file = paste0("Exponential Model dnds vs protist with exp priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 53){
  print("Exponential Model dnds vs protist with flat priors and with random IDs.Rdata")

  model.wpro.flat.rand<- MCMCglmm(w_ratios_Value ~ Protist,
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
  save(model.wpro.flat.rand, file = paste0("Exponential Model dnds vs protist with flat priors and with random IDs.Rdata"))
  print("Done")
}

if (k == 54){
  print("Exponential Model dnds vs protist with iw priors and with random IDs.Rdata")

  model.wpro.iw.rand<- MCMCglmm(w_ratios_Value ~ Protist,
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
  save(model.wpro.iw.rand, file = paste0("Exponential Model dnds vs protist with iw priors and with random IDs.Rdata"))
  print("Done")
}
