## Script for MCMCGLMM on selection parameters

### Preparation of the dataset

Run  Preparation_of_datasets_for_MCMCglmm.R

### MCMCglmm models

Using the dataset files produced and the calibrated tree, the models
were performed with 

dNdS for All OGs:
Run_dndsMCOlow_models.sh and MCMCglmm_MCOlowdnds.R

dNdS for single copy OGs:
Run_dnds_models.sh and MCMCglmm_SCOdnds.R

k for single copy OGs:
Run_relax_models.sh and MCMCglmm_SCOk.R

### Extracting results from models

The results from all models were extracted using

MCMCglmm_results_extraction.R and R_script_Foodcategory_Socialcategory.R
