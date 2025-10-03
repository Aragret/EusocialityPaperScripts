################################################################################
#
#       Title:    Extract the results from MCMCglmm models
#       Project:  Eusociality
#       Author:   Aumont
#       Year:     2024
#
################################################################################

#-------------------------------------------------------------------------------
#----------- Introduction ------------------------------------------------------
#-------------------------------------------------------------------------------

#----------- Aim ---------------------------------------------------------------

# We extract automatically in a table the results for models of which
# the explanatory variables are either covariates of binomial (but not
# categorical)

#----------- Library importation -----------------------------------------------

library(tidyverse)
library(ape)
library(geiger)
library(MCMCglmm)


#----------- Working directory -------------------------------------------------

setwd("C:/Users/caumont/Documents/WORK/Dino_McMahon/Work/Lab work and projects/10. Eusociality/MCMCGLMM/_Scripts")

#-------------------------------------------------------------------------------
#----------- Result exportation ------------------------------------------------
#-------------------------------------------------------------------------------

#----------- Selection of the models  ------------------------------------------
#Here we first select the models that we want to extract the information.
# The models have all been saved in a directory, and named according to the 
# analysis performed. To extract the result for all models testing dNdS values
# of MCO, we give "depvar" the value "MCOdNds"

depvar ="MCOdNdS" #dNdS #relax

# filenames list the names of the models selected
filenames = list.files("results",
                       pattern = paste0("*",depvar,"*"),
                       full.names = T)

#----------- Dataframe of the results (initialisation) -------------------------

# Here we create a dataframe with all the columns that need to be filled

N_row = length(filenames)

df = data.frame(list("Tested_variable" = rep("NA", N_row),
                     "Prior"=rep("NA",N_row),
                     "Effective_size"=rep("NA",N_row),
                     "Auto_correlation"=rep("NA",N_row),
                     "DIC"=rep("NA",N_row),
                     "Intercept_LowHPD"=rep("NA",N_row),
                     "Intercept"=rep("NA",N_row),
                     "Intercept_HighHPD"=rep("NA",N_row),
                     "Estimate_LowHPD"=rep("NA",N_row),
                     "Estimate"=rep("NA",N_row),
                     "Estimate_HighHPD"=rep("NA",N_row),
                     "Philo_LowHPD"=rep("NA",N_row),
                     "Philo_median"=rep("NA",N_row),
                     "Philo_mode"=rep("NA",N_row),
                     "Philo_HighHPD"=rep("NA",N_row),
                     "Ortho_LowHPD"=rep("NA",N_row),
                     "Ortho_median"=rep("NA",N_row),
                     "Ortho_mode"=rep("NA",N_row),
                     "Ortho_HighHPD"=rep("NA",N_row),
                     "random_id_LowHPD"=rep("NA",N_row),
                     "random_id_median"=rep("NA",N_row),
                     "random_id_mode"=rep("NA",N_row),
                     "random_id_HighHPD"=rep("NA",N_row),
                     "residuals_LowHPD"=rep("NA",N_row),
                     "residuals_median"=rep("NA",N_row),
                     "residuals_mode"=rep("NA",N_row),
                     "residuals_HighHPD"=rep("NA",N_row)
))

#----------- Dataframe of the results (filling) --------------------------------

# for each model we extract the information
j=0
for (modfile in filenames){
  print(modfile)
  j = j + 1
  
  # load the model
  mod = load(modfile)
  model = eval(as.name(mod))
  
  # extract the formula
  formula_mod = as.character(model$Fixed$formula)[3]
  
  # this indicate the type of random effect in the model
  random_mod = sum(model$Random$nfl)
  
  # extract prior name
  prior_mod = strsplit(strsplit(modfile, "with")[[1]][2], " ")[[1]][2]
  
  # extract effect size
  eff_size = effectiveSize(model$Sol[, 1:model$Fixed$nfl, 
                          drop = FALSE])[[1]]
  
  # extract maximale autocorrelation 
  autocorrelation = autocorr(model$VCV)
  max_autocor = max(autocorrelation[5,,])
  
  # extract DIC 
  DIC_mod = model$DIC
  
  # extract information on intercept (value and CI)
  intercept = summary(model)$solutions[1,1]
  intercept.low = summary(model)$solutions[1,2]
  intercept.high = summary(model)$solutions[1,3]
  
  # extract information on estimate (value and CI)
  estimate = summary(model)$solutions[2,1]
  estimate.low = summary(model)$solutions[2,2]
  estimate.high = summary(model)$solutions[2,3]
  
  # extract information on random_id (value and CI)
  if(random_mod == 3){
    random_id.VCV = model$VCV[,'random_id']
    }else{
    random_id.VCV = 0
    }
  random <- random_id.VCV/(phylo.VCV + ortho.VCV + random_id.VCV + residual.VCV)
  random.med = median(random)
  random.mode = as.numeric(posterior.mode(random))
  random.low = HPDinterval(random)[1]
  random.high = HPDinterval(random)[2]
  
  # extract information on phylogenetic effect (value and CI)
  phylo.VCV = model$VCV[,'Species']
  phylo <- phylo.VCV/(phylo.VCV + ortho.VCV + random_id.VCV + residual.VCV)
  phylo.med = median(phylo)
  phylo.mode = posterior.mode(phylo)
  phylo.low = HPDinterval(phylo)[1]
  phylo.high = HPDinterval(phylo)[2]
  
  # extract information on Orthogroup (value and CI)
  ortho.VCV = model$VCV[,'Orthogroup']
  ortho <- ortho.VCV/(phylo.VCV + ortho.VCV + random_id.VCV + residual.VCV)
  ortho.med = median(ortho)
  ortho.mode = as.numeric(posterior.mode(ortho))
  ortho.low = HPDinterval(ortho)[1]
  ortho.high = HPDinterval(ortho)[2]
  
  # extract information on residuals (value and CI)
  residual.VCV = model$VCV[,'units']
  residual <- residual.VCV/(phylo.VCV + ortho.VCV + random_id.VCV + residual.VCV)
  residual.med = median(residual)
  residual.mode = as.numeric(posterior.mode(residual))
  residual.low = HPDinterval(residual)[1]
  residual.high = HPDinterval(residual)[2]
  
  #fill up a raw with the information of the model
  df[j,] = c(formula_mod,
             prior_mod,
             eff_size,
             max_autocor,
             DIC_mod,
             intercept.low,
             intercept,
             intercept.high,
             estimate.low,
             estimate,
             estimate.high,
             phylo.low,
             phylo.med,
             phylo.mode,
             phylo.high,
             ortho.low,
             ortho.med,
             ortho.mode,
             ortho.high,
             random.low,
             random.med,
             random.mode,
             random.high,
             residual.low,
             residual.med,
             residual.mode,
             residual.high
             )
}
#check dataframe
str(df)

#----------- Dataframe of the results (saving) --------------------------------

write.table(df, paste0("Results_",depvar,".txt"), quote=F, dec=",", sep = "\t")

# For each type of analysis, the results have been copy pasted in the excel
# table Table_all_results.xlsx


