setwd('M:/Share/PhD/eusociality_paper/analyses/cpg_oe/')

library(dplyr)
library(stats)
library(tidyr)
library(ggplot2)

sc_hogs = read.table('SingleCopyHOGs', sep = '\t', header = FALSE) %>%
  select(c(V1, V4:V32)) %>%
  pivot_longer(!V1) %>%
  select(-name)

names(sc_hogs) = c('HOG', 'Name')

cpg_oe = read.table('all_CpGoe.txt', sep = '\t', header = TRUE) %>%
  filter(Name != 'Name') %>%
  mutate(Species = substr(Name, 1, 4)) %>%
  right_join(sc_hogs) %>%
  filter(!is.na(CpGoe))

cpg_oe$CpGoe = as.numeric(cpg_oe$CpGoe)

cpg_oe_wide = cpg_oe %>%
  # filter(!is.na(Species)) %>%
  select(c(Species, CpGoe, HOG)) %>%
  pivot_wider(names_from = HOG, values_from = CpGoe)


pca = prcomp(cpg_oe_wide[,-1])


explained_variance_ratio <- summary(pca)[["importance"]]['Proportion of Variance',]
explained_variance_ratio <- 100 * explained_variance_ratio
explained_variance_ratio

# PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11   PC12   PC13   PC14   PC15   PC16 
# 23.385 13.150  9.789  6.537  6.166  6.065  4.676  3.375  2.414  2.051  1.882  1.793  1.697  1.621  1.559  1.450 
# PC17   PC18   PC19   PC20   PC21   PC22   PC23   PC24   PC25   PC26   PC27   PC28   PC29 
# 1.400  1.366  1.271  1.213  1.146  1.127  1.035  0.966  0.900  0.768  0.727  0.472  0.000

components <- pca[["x"]]
components <- data.frame(components)
components = cbind(cpg_oe_wide$Species, components)

names(components)[1] = 'Species'

components_groups = components %>%
  mutate(
    Groups = case_when(
    Species %in% c('BGER', 'Bori', 'KAJ4', 'evm.') ~ 'Cockroaches',
    Species %in% c('Cmer', 'Cpun') ~ 'Cryptocercus',
    Species %in% c('Hsjo', 'Znev', 'Ncas', 'Kfla', 'PRsi', 'Cbre', 'Isch') ~ 'TermitesLinear',
    Species %in% c('Mdar', 'Dlon', 'Rfla', 'Hten', 'Ctes', 'Cges', 'Ofor', 'Mnat', 'Ssph', 'Aunk', 'Apac', 'Nluj', 'Csp4', 'Ntar', 'Pred', 'Munk') ~
      'TermitesBifurcated'
  ))

ggplot(components, aes(PC1, PC2, colour = as.factor(Species))) +
  geom_point()

ggplot(components_groups, aes(PC1, PC2, colour = as.factor(Groups))) +
  geom_point(size = 3) +
  theme_classic() +
  xlab('PC1 (23.38%)') +
  ylab('PC2 (13.15%)')

write.table(components[, c(1:3)], 'PCA.tsv', quote = F, sep = '\t', row.names = F)
