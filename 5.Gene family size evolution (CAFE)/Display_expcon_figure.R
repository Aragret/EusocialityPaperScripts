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
In this particular script, we aim to colour the branches from a 
phylogenetic tree with a gradient of colour following the number of 
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

# set to current wd

# ----- [1.1] load data -----------------------------------------------


full_data = read.csv2("../Data/dataset_for_figure.csv")


# -------- figures -------------------------------------------------------------

gradient_gain_loss = (c("cyan","cyan","cyan2","cyan2","cyan3","cyan4",
  'grey30'  ,
  "red4","red3",'red2','red2','red','red'))
legend_title_gainloss ="Average\ngene family\nexpansion\nand\ncontraction\nper Mya"
lim_gain_loss = c(-15,15)
ggsave_filename_gain_loss = "../Result/Figure4 CAFE tree gain-loss - events2.png"


contraction = full_data$contraction
expansion = full_data$expansion
gain_loss = (full_data$expansion - full_data$contraction) / (full_data$branch_size*100) # per Mya
gain_loss[30]=0
# remained = (0.5+(39-full_data$neutral)/39)
treeimg<-  ggtree(tree, aes(color=gain_loss), size=1) +
  scale_color_gradientn(colours = gradient_gain_loss,
                        name=legend_title_gainloss,
                        limits = lim_gain_loss,n.breaks = 2)+
  #scale_color_viridis(name='gain_loss', option="plasma") +
  geom_tiplab(fontface = 3,hjust = -0.15, size=3.4) +
  theme(legend.position=c(.26,.8), plot.margin = margin(6,120,6,6))+
  coord_cartesian(clip = "off")+
  #geom_point2(aes(subset=node==74), color='darkgreen', size=5)+
  geom_nodelab(aes(x=branch, 
                   label=paste0(expansion,"/",contraction,"/",round(gain_loss,2))),
               size = 1.6,
               vjust =#-0.5 
                       c(-.5, -.5, -.5, -.5, -.5, -.5, -.5, -.5, -.5, -.5,
                         -.5, -.5, -.5, -.5, -.5, -.5, 1.5, -.5, -.5, -.5,
                         1.5, -.5, -.5, 3.8, -.5, -.5, -.5, -.5)
               , 
               hjust = #1
                       c(1.2, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.8, 0.7, 0.6,
                         0.8, 0.8, 0.8, 0.6, 0.9, 1.0, 0.7, 0.4, 0.8, 0.7,
                         .95, 0.6, 1.0, 0.5, 0.9, 1.1, 1.0, 0.7)
               )+
  geom_tiplab(aes(x=branch, 
                  label=paste0(expansion,"/",contraction,"/",round(gain_loss,2)),colour=gain_loss),
              size = 1.6,
              #c(rep(alpha(1,1),30),"cyan2",rep(alpha(1,1),3),"cyan2",rep(alpha(1,1),12)),
              vjust = 
                c(-.5, -.5, -.5, -.5, -.5, -.5, -.5, -.5, -.5, -.5,
                  -.5, -.5, -.5, -.5, -.5, -.5, -.5, 1.5, -.5, -.5,
                  -.5, -.5, -.5, -.5, -.5, -.5, -.5, -.5, -.5),
              hjust = 0.5)+
  theme(legend.title = element_text(hjust = .5, size = 10.5), 
        legend.ticks = element_blank(), legend.title.position = "left")
treeimg


ggsave(file=ggsave_filename_gain_loss, plot=treeimg,
       width = 7, height = 7.5,dpi = 660)




