# Script to subset Melihs dN/dS data by caste-biased genes
# from:
# (/global/students/research/yldz/statisticalResults/selection/Codeml/median-dnds-tree.R).

library(ape)
# library(geiger)
# library(nlme)
library(phytools)
library(ggtree)
library(ggplot2)
# library(viridis)
# library(readr)
# library(cowplot)
library(dplyr)
library(tidyr)
library(reshape2)
library(lme4)

# setwd("~/PhD/selection")

tree <- ape::read.tree("230621_27610_99_unfilt_GTR_aln.nexus.treefile.tree")

codeml_tree = treeio::read.codeml_mlc("HOG0012341_free.out")@phylo

edge_codeml = data.frame(codeml_tree$edge)
colnames(edge_codeml) <- c("parent", "node")

# edge <- data.frame(tree$edge, edge_num = 1:nrow(tree$edge))
# colnames(edge) <- c("parent", "node", "edge_num")

dnds <- read.table("dnds.txt", sep = "\t", header = TRUE) %>%
  mutate(newCol = "dNdS") %>%
  pivot_wider(names_from = newCol, values_from = dNdS)

ds <- read.table("ds.txt", sep = "\t", header = TRUE) %>%
  mutate(newCol = "dS") %>%
  pivot_wider(names_from = newCol, values_from = dS)

dnds_matrix = t(matrix(unlist(dnds$dNdS), ncol = nrow(dnds)))
dnds_matrix[1, ] == dnds[1, ]$dNdS[[1]]

ds_matrix = t(matrix(unlist(ds$dS), ncol = nrow(ds)))

filter_dS = which(ds_matrix > 3)
ds_matrix[filter_dS]
dnds_matrix[filter_dS]

dnds_matrix[filter_dS] = NA
which(is.na(dnds_matrix)) == which(ds_matrix > 3)

dnds$median_dnds = apply(dnds_matrix, 1, function(x) median(x[x < 10], na.rm = TRUE))

edge_codeml = edge_codeml %>%
  mutate(Branch = paste(as.character(parent), as.character(node), sep = "..")) %>%
  full_join(dnds[, -2])

dnds_tree = full_join(codeml_tree, edge_codeml)

pdf("labelled_tree.pdf")
ggtree(dnds_tree, branch.length = "none") + geom_text(aes(label=node), hjust=-.3) + labs(title = "CodeML tree")
ggtree(dnds_tree, branch.length = "none") + geom_tiplab()
dev.off()


test = codeml_tree
test$edge = edge_codeml

p = ggtree(test, aes(color = median_dnds), branch.length = "none", size =3) +
  geom_tiplab(inherit.aes=T) + ggtitle("Median dN/dS across branches") +
  geom_text(aes(label=label), size=3, hjust = -1)



###GLMM
library(lme4)
library(reshape2)
colnames(dnds) <- c("mmus","hgla","fdam","X4","cpor","clan","odeg","X8","X9","X10")
dnds$OG <- rownames(dnds)
dnds.m <- melt(dnds[,c(1:3,5:7,11)])
colnames(dnds.m)[2:3] <- c("species","dnds")
sociality <- data.frame(species = c("mmus","hgla","fdam","cpor","clan","odeg"),sociality = c(0,2,1,0,0,0))
#sociality <- data.frame(species = c("mmus","hgla","fdam","cpor","clan","odeg"),sociality = c("N","ES","ES",rep("N",3)))
age <- data.frame(species = c("mmus","hgla","fdam","cpor","clan","odeg"), age = c(1,30,20,6,10,7))
colony <- data.frame(species = c("mmus","hgla","fdam","cpor","clan","odeg"), colony = c(3,295,40,10,100,10))
dnds.m <- merge(dnds.m,sociality,by="species")
dnds.m <- merge(dnds.m, age,by="species")
dnds.m <- merge(dnds.m, colony,by="species")
dnds.glmm <- glmer(dnds~ colony + age + sociality + (1|OG) + (1|species), family=Gamma(link = "log"), data=dnds.m[dnds.m$dnds<10,])
summary(dnds.glmm)


### boxplots


dnds_matrix[which(dnds_matrix > 10)] = NA

df_matrix = as.data.frame(dnds_matrix)
df_matrix$Branch = dnds$Branch

df = melt(df_matrix, id.vars = "Branch")

# get all descendants of the node where bifurcated dev emerged

bifdev = data.frame(node = phytools::getDescendants(codeml_tree, 56))

phytools::getParent(codeml_tree, c(57))

bifdev$parent = lapply(bifdev$node, function(x) getParent(codeml_tree, c(x)))

bifdev = bifdev %>% 
  mutate(Branch = paste(parent, node, sep = "..")) %>%
  select(Branch)

rest = data.frame(Branch = c("50..28", "83..15"))

bifdev = rbind(bifdev, rest) %>%
  mutate(BifurDev = 1)

pdf("boxplots.pdf", width = 50)
ggplot(df, aes(Branch, value)) +
  geom_boxplot() + ylim(0, 1)
dev.off()

df = full_join(df, bifdev)

df$BifurDev[is.na(df$BifurDev)] = 0

df$BifurDev = as.factor(df$BifurDev)

names(df) = c('Branch', 'HOG', 'dnds', 'BifurDev')

dnds.glmm <- glmer(dnds ~ BifurDev + (1|HOG), family=Gamma(link = identity), df)

summary(dnds.glmm)

df = df %>%
  mutate(Sociality = as.factor(case_when(
    Branch == "48..6" ~ "Solitary",
    Branch %in% c("49..10", "48..49") ~ "Subsocial",
    BifurDev == 1 ~ "BifurDev",
    .default = "LinearDev"
  )))

dnds.glmm <- glmer(dnds ~ Sociality + (1|HOG), family=Gamma(link = identity), df[!is.na(df$dnds),])
summary(dnds.glmm)

summary(lm(dnds ~ Sociality, df))

summary(lm(dnds ~ relevel(Sociality, ref="Solitary"), df))

summary(lm(dnds ~ relevel(Sociality, ref="Subsocial"), df))

summary(lm(dnds ~ relevel(Sociality, ref="LinearDev"), df))
