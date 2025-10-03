################################################################################
#
#             Title: Colour expanded/contracted branches from CAFE 
#           Project: Eusociality
#           Authors: C. Aumont
#              Year: 2024
#
################################################################################


# ------------------------------------------------------------------------------
# ----- [0] Description --------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- [0.1] Purpose ----------------------------------------------------------

'
In this particular script, we aim to prepare the dataset to colour branches 
from a phylogenetic tree with a gradient of colour following the number of 
expansion/contraction per branch length
'


# ----- [0.2] Dependancies: ----------------------------------------------------

'The following libraries are required: ape, phytools, ggplot2, dplyr,
ggtree, viridis'
'The following input files are:
tcal_blattodea_tree.tre
selected_model_for_small_hog/Gamma_change.tab
selected_model_for_large_hog/Gamma_change.tab
selected_model_for_small_hog/Gamma_branch_probabilities.tab
selected_model_for_large_hog/Gamma_branch_probabilities.tab
'

# ------------------------------------------------------------------------------
# ------ [1] Introduction ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- [1.1] Load packages ----------------------------------------------------
library(ape)
library(phytools)
library(dplyr)
library("ggtree")
library(ggplot2)
library(viridis)


# ----- [1.2] set work directory -----------------------------------------------

# set to current workdir

# ----- [1.1] load tree and extract ages -----------------------------------------------


tree=read.tree("../Data/tcal_blattodea_tree.tre")
plot(tree)
nodelabels()
indice = sapply(1:57,function(x,y) which (y==x),
           y=tree$edge[,2])

data1=setNames(tree$edge.length[unlist(indice)],c(tree$tip.label,as.character(c(31:57))))

data=c(data1[1:29],"30"=0,data1[30:56])

df_branchsize=data.frame("branch_size" = data)

write.table(df_branchsize, file = "../Result/Ages_branchsize_eusociality.txt",quote = F,
           dec = ",", sep = "\t",row.names = T)


# ----- [1.1] con and exp count -----------------------------------------------

#---------------------------- change table creation ----------------------------

change.tabS0 = read.table("../Data/selected_model_for_small_hog/Gamma_change.tab",
                          header = T,dec = ".",na.strings = "<NA>")
rownames(change.tabS0)=change.tabS0$FamilyID
change.tabS=change.tabS0[,-1]
head(change.tabS)
change.tabL0 = read.table("../Data/selected_model_for_large_hog/Gamma_change.tab",
                          header = T,dec = ".",na.strings = "<NA>")
rownames(change.tabL0)=change.tabL0$FamilyID
change.tabL=change.tabL0[,-1]
head(change.tabL)
change.tab = rbind.data.frame(change.tabS, change.tabL)
tail(change.tab,n = 15)


#--------------------- Evolution probability table creation --------------------


proba.tabS0 = read.table("../Data/selected_model_for_small_hog/Gamma_branch_probabilities.tab",
                         header = T,dec = ".",na.strings = "<NA>")
rownames(proba.tabS0)=proba.tabS0$FamilyID
proba.tabS=proba.tabS0[,-1]
head(proba.tabS)
proba.tabL0 = read.table("../Data/selected_model_for_large_hog/Gamma_branch_probabilities.tab",
                         header = T,dec = ".",na.strings = "<NA>")
rownames(proba.tabL0)=proba.tabL0$FamilyID
proba.tabL=proba.tabL0[,-1]
head(proba.tabL)
proba.tab = rbind.data.frame(proba.tabS, proba.tabL)
tail(proba.tab,n = 15)

proba = t.data.frame(proba.tab)
proba = data.frame(proba)
change = t.data.frame(change.tab)
change = data.frame(change)
change.tab.sig = subset.data.frame(change, select = colnames(proba))
cc = t.data.frame(change.tab.sig)
tail(cc)
dim(cc)
change.tab = cc

n_node = as.numeric(dim(change.tab)[2])
#head(pre_data)
data.change=data.frame(row.names = row.names(df))
data.change$neutral = rep(0,n_node)
data.change$expansion = rep(0,n_node)
data.change$contraction = rep(0,n_node)


for (i in (1:dim(change.tab)[2])){
  for (u in (1:dim(change.tab)[1])){
    if (proba.tab[u,i] > 0.05){change.tab[u,i]=0}
  }
  for (j in change.tab[,i]){
    
    if (j==0){
      data.change$neutral[i]=data.change$neutral[i]+1
    }
    if (j>0){
      data.change$expansion[i]=data.change$expansion[i]+1

    }
    if (j<0){
      data.change$contraction[i]=data.change$contraction[i]+1

    }
  }
}


write.csv2(change.tab,"../Result/table_of_expansioncontraction_per_hog_on_tree.csv",
           quote = F,row.names = T)

write.csv2(subset(change.tab,row.names(change.tab) %in% c("N0.HOG0000380",
                                                          "N0.HOG0001851",
                                                          "N0.HOG0002046",
                                                          "N0.HOG0002649",
                                                          "N0.HOG0004986",
                                                          "N0.HOG0000937",
                                                          "N0.HOG0001513",
                                                          "N0.HOG0002943",
                                                          "N0.HOG0006820",
                                                          "N0.HOG0000758"
                                                          )),
           "../Result/table_of_expansioncontraction_per_hog_on_tree_subset.csv",
           quote = F,row.names = T)

proba.tab2 =proba.tab
for (i in (1:dim(change.tab)[2])){
  for (u in (1:dim(change.tab)[1])){
    if (proba.tab[u,i] > 0.05){proba.tab2[u,i]=1}
  }
}
write.csv2(subset(proba.tab2,row.names(proba.tab) %in% c("N0.HOG0000380",
                                                          "N0.HOG0001851",
                                                          "N0.HOG0002046",
                                                          "N0.HOG0002649",
                                                          "N0.HOG0004986",
                                                          "N0.HOG0000937",
                                                          "N0.HOG0001513",
                                                          "N0.HOG0002943",
                                                          "N0.HOG0006820",
                                                          "N0.HOG0000758"
)),
"../Result/table_of_pval_expcont_per_hog_on_tree_subset.csv",
quote = F,row.names = T)

### significant gene, expansion or contraction

test = change.tab
cont_list=c()
expa_list=c()
evol_list=c()
for (i in (1:dim(change.tab)[1])){
  cont=0
  expa=0
  for (j in (1:dim(change.tab)[2])){
    if (test[i,j] != 0){
    test[i,j] = change.tab[i,j]/abs(change.tab[i,j])
    if (test[i,j] < 0){cont = cont +1}
    if (test[i,j] > 0){expa = expa +1}
    }
  }
  cont_list = append(cont_list, cont)
  expa_list = append(expa_list, expa)
  evol_list = append(evol_list, expa - cont)
}
evol_df = data.frame(expansions = expa_list,
                   contractions=cont_list,
                   evolution = evol_list)
resexp = 0
rescon = 0
for (i in evol_df$evolution){
  if (i > 0) {resexp = resexp + 1}
  if (i < 0) {rescon = rescon + 1}
}
resexp
rescon
###

data.change$total=data.change$neutral+data.change$expansion+data.change$contraction
head(data.change)
tail(data.change)

# Merge datasets 

full_data = cbind.data.frame(df_branchsize,data.change)

write.csv2(full_data, "../Result/dataset_for_figure.csv")






