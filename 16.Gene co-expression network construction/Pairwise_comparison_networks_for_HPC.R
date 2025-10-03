# For HPC

# ----- [1] Setting up ---------------------------------------------------------

# ----- [1.1] uploading libraries ----------------------------------------------

print("uploading libraries...")
library(CoDiNA, lib.loc = '/home/cedria95/ca3321fu/R/x86_64-pc-linux-gnu-library/4.2')

# ----- [1.2] Setting working directory ----------------------------------------

print("Setting working directory...")
setwd("/scratch/cedria95/gene_network_analyses/data_preparation/4_wTO")

# ----- [1.3] Uploading arguments ----------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
print("Working with species")
spe = as.character(args)
print(spe)

# ----- [2] Networks -----------------------------------------------------------

# ----- [2.1] Uploading networks -----------------------------------------------

print("loading Networks...")
Networkfilelist = list.files("result_casteNormalised",      # to change if use new ORs annotation
                         pattern = paste0(spe,"_Wo|",spe,"_Nm|",spe,"_Er|",spe,"_Pr|",spe,"_Al|",spe,"_Ne"),
                         full.names = T)

Network_list = list()
network_codes = c()
# Networkfile=Networklist[1]
for (Networkfile in Networkfilelist){
  print(Networkfile)
  caste = unlist(strsplit(x = unlist(strsplit(Networkfile, "_"))[5],".R"))[1]
  Networkname = load(Networkfile)
  Network = eval(as.name(Networkname))
  Network$Node.1=as.character(Network$Node.1)
  Network$Node.2=as.character(Network$Node.2)
  Network_spe_caste = paste0("Network_",spe,"_",caste)
  network_codes = append(network_codes,Network_spe_caste)
  assign(Network_spe_caste, Network)
  Network_list=append(Network_list,list(eval(as.name(Network_spe_caste))))
}
print("done")


# ----- [2.2] Comparing Networks -----------------------------------------------

print("Comparing networks...")
Diff_Caste = MakeDiffNet(Data = Network_list,
                         Code = network_codes,
                         cutoff = 0.33,
                         stretch = TRUE)
print("done")


# ----- [2.3] Clustering genes -------------------------------------------------

print("Clustering genes...")
int_C = quantile(Diff_Caste$Score_internal, 0.5)
ext_C = quantile(Diff_Caste$Score_Phi, 0.5)

Nodes_Groups = ClusterNodes(DiffNet = Diff_Caste, 
                            cutoff.external = ext_C, 
                            cutoff.internal = int_C)

Group_table = table(Nodes_Groups$Phi_tilde)


# ----- [2.4] saving clustering ------------------------------------------------

print("saving Group table...")
csv_name=paste0("Codina_pairwise_WORE_beforeBC/PairwiseWoReBeforeBC_Group_table_normalised_",spe,".csv")
write.csv2(x = Group_table , file =  csv_name)
print(paste0(csv_name,"saved"))

print("saving Node groups...")
csv_name=paste0("Codina_pairwise_WORE_beforeBC/PairwiseWoReBeforeBC_Node_groups_normalised_",spe,".csv")
write.csv2(x = Nodes_Groups , file =  csv_name)
print(paste0(csv_name,"saved"))
