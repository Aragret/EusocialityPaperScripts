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
# the explanatory variables are categorical

#----------- Library importation -----------------------------------------------

library(tidyverse)
library(ape)
library(geiger)
library(MCMCglmm)


# ---------- Set work directory ------------------------------------------------

setwd("C:/Users/caumont/Documents/WORK/Dino_McMahon/Work/Lab work and projects/10. Eusociality/MCMCGLMM/_Scripts")

#-------------------------------------------------------------------------------
#----------- Load datasets -----------------------------------------------------
#-------------------------------------------------------------------------------
# Each model is loaded and the summary is saved automatically

# soc or food
filenames = list.files("results", pattern = "relax", full.names = T) 

for (modfile in filenames){
  print(modfile)
  mod = load(modfile)
  model = eval(as.name(mod))
  print(summary(model))
  
  # model results are saved in a csv file
  write.csv2(summary(model)$solution, file = paste0("Other_soc_",modfile,".csv"))
}

# For each type of analysis, the representative results have been copy pasted 
# in the excel table Table_all_results.xlsx
