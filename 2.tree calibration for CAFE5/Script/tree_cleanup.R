################################################################################
#
#             Title: Blattodea tree preparation for CAFE
#           Project: Eusociality
#           Authors: C. Aumont
#              Year: 2025
#
################################################################################


# ------------------------------------------------------------------------------
# ----- [0] Description --------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- [0.1] Purpose ----------------------------------------------------------

'
This script is used to modify the tree from the Eusociality
project so that it can be used in CAFE5 for gene family size evolution.

In this script, the original tree from the Eusociality project is
named "mcmctree1_wag_runA_1Msamples_figtree.tre". The output tree is called
"tcal_blattodea_tree.tre" and used in further analysis with CAFE5.

We import the tree and check for binary, root presence, and ultrametry. 
Ultrametry and tip names are adjusted. The resulting tree is saved.
'


# ----- [0.2] Dependancies and input files -------------------------------------

'The following libraries are required: ape, geiger, phytools'

'The following files are required:
mcmctree1_wag_runA_1Msamples_figtree.tre
Blattodea_names.csv
'

# ------------------------------------------------------------------------------
# ------ [1] Introduction ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- [1.1] Load packages ----------------------------------------------------

library(ape)
library(geiger)
library(phytools)


# ----- [1.2] Set directory ----------------------------------------------------

path_cafe = "C:/Users/caumont/Documents/WORK/Dino_McMahon/Work/Lab work and projects/10. Eusociality/CAFE/"
setwd(paste0(path_cafe,"1. tree"))

# ----- [1.3] Load tree --------------------------------------------------------

# Read in the tree
tree1 <- read.tree("mcmctree1_wag_runA_1Msamples_figtree.tre")

# Look at the tree summary
tree1

# ----- [1.4] Display the tree -------------------------------------------------

# Plot the tree with small labels
par(mfrow=c(1,1))
plot(tree1, cex = 0.6, typ = "phylogram", no.margin = TRUE)

# ----- [1.5] Check tree properties --------------------------------------------

# Check if tree parameters are correct
is.binary(tree1) # True
is.rooted(tree1) # True
is.ultrametric(tree1) # False

# ------------------------------------------------------------------------------
# ----- [2] Ultrametry ---------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- [2.1] Force ultrametry -------------------------------------------------

# We use force.ultrametric from phytools with extend method to round up the edge
Utree = force.ultrametric(tree1, method = "extend")
plot(Utree)


# ----- [2.3] Check changes ----------------------------------------------------

# We can plot the tree and compare with the previous plot.
plot(Utree, cex = 0.6, typ = "phylogram", no.margin = TRUE)
plot(tree1, cex = 0.6, typ = "phylogram", no.margin = TRUE)

# More accurately, we can compare the edges before and after forcing ultrametry 
# to be sure that it was indeed a problem of rounding the branch length
unclass(Utree)$edge.length - unclass(tree1)$edge.length

# 
# > unclass(Utree)$edge.length - unclass(tree1)$edge.length
# [1] 0e+00 3e-06 3e-06 0e+00 0e+00 3e-06 3e-06 0e+00 0e+00 2e-06 2e-06 0e+00 2e-06 0e+00
# [15] 0e+00 1e-06 1e-06 0e+00 0e+00 2e-06 0e+00 2e-06 0e+00 2e-06 2e-06 0e+00 1e-06 0e+00
# [29] 1e-06 0e+00 0e+00 1e-06 0e+00 1e-06 0e+00 1e-06 1e-06 0e+00 0e+00 0e+00 0e+00 0e+00
# [43] 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 0e+00 1e-06 1e-06


# ------------------------------------------------------------------------------
# ----- [3] Tips names ---------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- [3.1] Load name dataset ------------------------------------------------


Blattodea_names <- read.csv("Blattodea_names.csv",
                           header = T,
                           sep = ";",
                           dec = ",",
                           stringsAsFactors = T)
summary(Blattodea_names)


# ----- [3.1] Check identity ---------------------------------------------------

# Check that termite names are identical as in the tree 
check <- name.check(phy = Utree,
                    data = Blattodea_names,
                    data.names = Blattodea_names$Species)
check


# ----- [3.2] Correct identity -------------------------------------------------

Utree$tip.label <- gsub("PRsi", "Psim", Utree$tip.label)
Utree$tip.label <- gsub("Csp4", "Cunk", Utree$tip.label)

# ------------------------------------------------------------------------------
# ------- [4] Save tree --------------------------------------------------------
# ------------------------------------------------------------------------------

write.tree(Utree, file = "tcal_blattodea_tree.tre")



