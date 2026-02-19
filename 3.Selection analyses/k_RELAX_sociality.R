library(ggplot2)
library(ape)
library(dplyr)
library(cowplot)
library(nlme)
library(lme4)
library(geiger)
library(phytools)
library(gghighlight)
library(ggtree)
library(aplot)
library(stringr)
library(gridExtra)


# setwd('/global/scratch2/amikhail/termite_genomes/eusociality_paper/relax/')
# setwd("~/EusocialityPaper/relax/")
setwd('M:/Share/eusociality_paper/selection')

df = read.csv('relax/singe_branch_relax.csv') %>%
  group_by(HOG) %>%
  mutate(p.adjHOG = p.adjust(p.value, method = 'BH')) %>%
  ungroup() %>%
  group_by(Branch) %>%
  mutate(p.adjBranch = p.adjust(p.value, method = 'BH'))

df$p.adj = p.adjust(df$p.value, method = 'BH')

######################

df$Species = ''

df[df$Branch == 'branch_1',]$Species = 'Dpun'
df[df$Branch == 'branch_2',]$Species = 'Bger'
df[df$Branch == 'branch_4',]$Species = 'Bori'
df[df$Branch == 'branch_5',]$Species = 'Pame'
df[df$Branch == 'branch_7',]$Species = 'Cmer'
df[df$Branch == 'branch_8',]$Species = 'Cpun'
df[df$Branch == 'branch_10',]$Species = 'Mdar'
df[df$Branch == 'branch_11',]$Species = 'Znev'
df[df$Branch == 'branch_12',]$Species = 'Hsjo'
df[df$Branch == 'branch_14',]$Species = 'Kfla'
df[df$Branch == 'branch_15',]$Species = 'Ncas'
df[df$Branch == 'branch_16',]$Species = 'Cbre'
df[df$Branch == 'branch_17',]$Species = 'Isch'
df[df$Branch == 'branch_21',]$Species = 'Dlon'
df[df$Branch == 'branch_22',]$Species = 'PRsi'
df[df$Branch == 'branch_23',]$Species = 'Rfla'
df[df$Branch == 'branch_24',]$Species = 'Hten'
df[df$Branch == 'branch_25',]$Species = 'Ctes'
df[df$Branch == 'branch_26',]$Species = 'Cges'
df[df$Branch == 'branch_30',]$Species = 'Ssph'
df[df$Branch == 'branch_31',]$Species = 'Mnat'
df[df$Branch == 'branch_32',]$Species = 'Ofor'
df[df$Branch == 'branch_35',]$Species = 'Apac'
df[df$Branch == 'branch_36',]$Species = 'Aunk'
df[df$Branch == 'branch_38',]$Species = 'Munk'
df[df$Branch == 'branch_39',]$Species = 'Pred'
df[df$Branch == 'branch_40',]$Species = 'Ntar'
df[df$Branch == 'branch_41',]$Species = 'Csp4'
df[df$Branch == 'branch_42',]$Species = 'Nluj'

Termites_linear = c("Cbre", "Isch", "Kfla", "Ncas", "PRsi", "Znev", "Hsjo")

Termites_bifur = c("Apac", "Aunk", "Cges", "Csp4", "Ctes", "Dlon", "Hten", "Mdar", "Mnat", "Munk",
                   "Nluj", "Ntar", "Ofor", "Pred", "Rfla", "Ssph")

Cryptocercus = c("Cmer", "Cpun")

Roaches = c("Bger", "Bori", "Dpun", "Pame")

##################

### adding aBSREL

absrel = read.table("absrel/all_significant") %>%
  mutate(Species = sub(",", "", V2)) %>%
  filter(!str_detect(Species, "Node"))

absrel$p.adj = p.adjust(absrel$V5, method = 'BH')

absrel_pcorr = filter(absrel, p.adj < 0.05)
# all significant

absrel_counts = as.data.frame(table(absrel$Species))
names(absrel_counts) = c("Species", "Positive")

### adding dNdS ratios

dnds = read.csv("codeml/dNdS_allHOGs.csv") %>%
  select(dn_Value, ds_Value, w_ratios_Value, Species.y, Orthogroup) %>%
  mutate(
    Species = case_when(
      Species.y == "BGER" ~ "Bger",
      Species.y == "evm." ~ "Dpun",
      Species.y == "KAJ4" ~ "Pame",
      TRUE ~ Species.y
    ),
    Social = factor(case_when(
      Species %in% Termites_linear ~ 'Termites_linear',
      Species %in% Termites_bifur ~ "Termites_bifur",
      Species %in% Cryptocercus ~ "Cryptocercus",
      Species %in% Roaches ~ "Solitary_cockroaches"
    ), levels = c('Solitary_cockroaches', 'Cryptocercus', 'Termites_linear', 'Termites_bifur'))
  )

dnds_sc = read.csv("codeml/GLMMtable_singlecopy.csv") %>%
  mutate(Social = factor(case_when(
    Species %in% Termites_linear ~ 'Termites_linear',
    Species %in% Termites_bifur ~ "Termites_bifur",
    Species %in% Cryptocercus ~ "Cryptocercus",
    Species %in% Roaches ~ "Solitary_cockroaches"
  ), levels = c('Solitary_cockroaches', 'Cryptocercus', 'Termites_linear', 'Termites_bifur'))
  ) %>%
  filter(ds_Value < 3 & dn_Value < 10)

###################################################################################

df_terminal_br = df[df$Species != '',] %>%
  mutate(Social = factor(case_when(
      Species %in% Termites_linear ~ 'Termites_linear',
      Species %in% Termites_bifur ~ "Termites_bifur",
      Species %in% Cryptocercus ~ "Cryptocercus",
      Species %in% Roaches ~ "Solitary_cockroaches"
    ), levels = c('Solitary_cockroaches', 'Cryptocercus', 'Termites_linear', 'Termites_bifur')))

table(df_terminal_br$Species)

ggplot(df_terminal_br, aes(Species, log(k))) +
  geom_violin() + ylim(c(-10, 10))

# write.csv(df_terminal_br, 'singe_branch_relax_species.csv', quote = F, row.names = F)

df_terminal_br_medians = df_terminal_br %>%
  group_by(Species) %>%
  summarise(k_median = median(k))

####### dNdS stats

by(dnds_sc$w_ratios_Value, dnds_sc$Social, summary)

# dnds_sc$Social: Solitary_cockroaches
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00010 0.02846 0.06537 0.11849 0.13037 9.75022 
# ------------------------------------------------------ 
#   dnds_sc$Social: Cryptocercus
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00010 0.06404 0.13638 0.20044 0.25719 3.29760 
# ------------------------------------------------------ 
#   dnds_sc$Social: Termites_linear
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0001  0.0741  0.1492  0.1879  0.2556  3.1075 
# ------------------------------------------------------ 
#   dnds_sc$Social: Termites_bifur
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00010 0.07873 0.17282 0.22822 0.30831 9.98742 

by(dnds$w_ratios_Value, dnds$Social, summary)

# dnds$Social: Solitary_cockroaches
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00010 0.03059 0.07739 0.18799 0.17050 9.94302 
# ------------------------------------------------------ 
#   dnds$Social: Cryptocercus
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00010 0.07475 0.17421 0.27873 0.33766 9.96421 
# ------------------------------------------------------ 
#   dnds$Social: Termites_linear
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00010 0.08258 0.18034 0.25927 0.32394 9.97543 
# ------------------------------------------------------ 
#   dnds$Social: Termites_bifur
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00010 0.08737 0.21161 0.31980 0.39998 9.92919 

####### Correlations

# cor.test(df_terminal_br_medians$Intensified, df_terminal_br_medians$Positive, method = "spearman")
# # rho 0.5942625, p-value 0.0006759
# 
# cor.test(df_terminal_br_medians$Intensified, df_terminal_br_medians$Relaxed, method = "spearman")
# 
# cor.test(df_terminal_br_medians$Positive, df_terminal_br_medians$Relaxed, method = "spearman")


####### Descr stats for the groups

relaxed_sign = df_terminal_br[df_terminal_br$p.value < 0.05 & df_terminal_br$k < 1,]
relaxed_sign_corr = df_terminal_br[df_terminal_br$p.adj < 0.05 & df_terminal_br$k < 1,]

relaxed_counts = as.data.frame(table(relaxed_sign$Species))
names(relaxed_counts) = c('Species', 'Relaxed')

relaxed_counts_corr = as.data.frame(table(relaxed_sign_corr$Species))
names(relaxed_counts_corr) = c('Species', 'Relaxed')

intense_sign = df_terminal_br[df_terminal_br$p.value < 0.05 & df_terminal_br$k > 1,]
intense_sign_corr = df_terminal_br[df_terminal_br$p.adj < 0.05 & df_terminal_br$k > 1,]

intense_counts = as.data.frame(table(intense_sign$Species))
names(intense_counts) = c('Species', 'Intensified')

intense_counts_corr = as.data.frame(table(intense_sign_corr$Species))
names(intense_counts_corr) = c('Species', 'Intensified')

df_terminal_br_medians = df_terminal_br_medians %>%
  mutate(
    Social = factor(case_when(
      Species %in% Termites_linear ~ 'Termites_linear',
      Species %in% Termites_bifur ~ "Termites_bifur",
      Species %in% Cryptocercus ~ "Cryptocercus",
      Species %in% Roaches ~ "Solitary_cockroaches"
    ), levels = c('Solitary_cockroaches', 'Cryptocercus', 'Termites_linear', 'Termites_bifur'))
  ) %>%
  left_join(intense_counts_corr) %>%
  left_join(relaxed_counts_corr) %>%
  left_join(absrel_counts)

pdf("../SelectionVis.pdf")

ggplot(df_terminal_br_medians, aes(Social, k_median)) +
  # geom_violin() +
  geom_jitter(width = 0.1)
  # gghighlight(Species %in% c("Mdar", "Dlon"))

ggplot(df_terminal_br_medians, aes(Social, Relaxed)) +
  # geom_violin() +
  geom_jitter(width = 0.1)
  # gghighlight(Species %in% c("Mdar", "Dlon"))

ggplot(df_terminal_br_medians, aes(Social, Intensified)) +
  # geom_violin() +
  geom_jitter(width = 0.1)
  # gghighlight(Species %in% c("Mdar", "Dlon"))

ggplot(df_terminal_br_medians, aes(Social, Positive)) +
  # geom_violin() +
  geom_jitter(width = 0.1) 
  # gghighlight(Species %in% c("Mdar", "Dlon"))


plot1 = ggplot(dnds, aes(Social, w_ratios_Value)) +
  geom_boxplot() +
  labs(title = "All OGs")

plot2 = ggplot(dnds_sc, aes(Social, w_ratios_Value)) +
  geom_boxplot() +
  labs(title = "Single-copy OGs")

plot3 = ggplot(dnds, aes(Social, ds_Value)) +
  geom_boxplot()

plot4 = ggplot(dnds_sc, aes(Social, ds_Value)) +
  geom_boxplot()

plot5 = ggplot(dnds, aes(Social, dn_Value)) +
  geom_boxplot()

plot6 = ggplot(dnds_sc, aes(Social, dn_Value)) +
  geom_boxplot()

grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 2)

plot7 = ggplot(dnds, aes(Social, w_ratios_Value)) +
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(0, 1))  +
  labs(title = "All OGs")

plot8 = ggplot(dnds_sc, aes(Social, w_ratios_Value)) +
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(0, 1))  +
  labs(title = "Single-copy OGs")

plot9 = ggplot(dnds, aes(Social, dn_Value)) +
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(0, 0.3))

plot10 = ggplot(dnds_sc, aes(Social, dn_Value)) +
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(0, 0.3))

plot11 = ggplot(dnds, aes(Social, ds_Value)) +
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(0, 3))

plot12 = ggplot(dnds_sc, aes(Social, ds_Value)) +
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(0, 2))

grid.arrange(plot7, plot8, plot9, plot10, plot11, plot12, ncol = 2)

dev.off()

by(df_terminal_br_medians, df_terminal_br_medians$Social, summary)

write.csv(df_terminal_br_medians, "RELAX_medians.csv")

###################################################################################

# tree <- ape::read.tree("Tree_NodeIDs.txt")
tree <- ape::read.tree("SpeciesTree_mod.txt")


tree$tip.label

ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

p <- ggtree(tree) %<+% df_terminal_br_medians +
  # coord_fixed(ratio = 1.1) + 
  geom_tiplab(align=TRUE)

p_rotated = rotate(p, 30)

# final = rotate(p, 30) + geom_tiplab(offset = 0.2) +
#   geom_tippoint(aes(color = log2(k_median)), size = 5) + 
#   # geom_nodepoint(aes(color = log2(k_median))) + 
#   guides(fill = "none") + # scale_size_continuous(range = c(1, 5)) + 
#   scale_color_continuous(type = "viridis") 

# p + geom_tiplab(offset = 0.2) +
  # geom_facet(panel = "Relaxed gene count", data = df_terminal_br_medians, geom = geom_col, 
  #            aes(x = Relaxed), orientation = 'y', width = .6)

dnds_plot = ggplot(dnds, aes(Species, w_ratios_Value, fill = Social)) + 
  geom_boxplot(outlier.shape = NA) + 
  # geom_hline(yintercept = 1, linetype="dashed", color = "red") +
  coord_flip(ylim = c(0, 1.2)) + theme_tree2() + 
  guides(fill = "none")

k_plot = ggplot(df_terminal_br, aes(Species, k, fill = Social)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_hline(yintercept = 1, linetype="dashed", color = "red") +
  coord_flip(ylim = c(0, 11)) + theme_tree2() + 
  guides(fill = "none")
  
relax_plot = ggplot(df_terminal_br_medians, aes(Species, Relaxed, fill = Social)) + 
  geom_col() + 
  coord_flip() + theme_tree2() +
  guides(fill = "none")

intens_plot = ggplot(df_terminal_br_medians, aes(Species, Intensified, fill = Social)) + 
  geom_col() + 
  coord_flip() + theme_tree2() +
  guides(fill = "none")

absrel_plot = ggplot(df_terminal_br_medians, aes(Species, Positive, fill = Social)) + 
  geom_col() + 
  coord_flip() + theme_tree2() +
  guides(fill = "none")

multi_panel = dnds_plot %>% insert_left(p, width = 3) %>% insert_right(k_plot) %>% insert_right(relax_plot) %>% insert_right(intens_plot) %>%
  insert_right(absrel_plot)

save_plot("SelectionTree.pdf", multi_panel,
          base_height = 5, base_width = 14)

# save_plot("~/EusocialityPaper/Figure_tree/SelectionTree.pdf", multi_panel,
#           base_height = 4, base_width = 10)

# save_plot('/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/figures/k_medians_2.pdf', final,
#           base_height = 4, base_width = 11)

###################################################################

### RELAX TW branches

relax_tw = read.table('relax/relax_table.tsv', header = TRUE, sep = '\t')

relax_tw_adj = relax_tw %>%
  mutate(p.adj = p.adjust(p.value, method = 'BH'))

relax_tw_relaxed = relax_tw_adj %>%
  filter(k < 1 & p.adj < 0.1)

relax_tw_int = relax_tw_adj %>%
  filter(k > 1 & p.adj < 0.1)

write.table(relax_tw_adj, 'relax_table_padj.tsv', sep = '\t', quote = F, row.names = F)

write.table(relax_tw_relaxed, 'relax_table_padj_relaxed.tsv', sep = '\t', quote = F, row.names = F)

write.table(relax_tw_int, 'relax_table_padj_int.tsv', sep = '\t', quote = F, row.names = F)


####################################################################################

### Correlating with the ontogeny metric

sociality_matrix = read.csv('/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/glmm/Copy_of_Matrix_of_ecological traits_13.06.2024 - Formated_matrix.csv', dec = ',') %>%
  select(Ontogeny, Nesting_type, Colony_size_log10_lower_range, Sociality_score, Sp_code, Sp_code_2) %>%
  filter(Sp_code_2 != '')

sociality_matrix$Species = sociality_matrix$Sp_code_2

df_sociality_matrix = 
  left_join(df_terminal_br, sociality_matrix)

df_sociality_matrix$Sociality_score = as.numeric(df_sociality_matrix$Sociality_score)

ggplot(df_sociality_matrix, aes(log(k))) +
  geom_density()

summary(glmer(log(k) ~ Sociality_score + (1|HOG) + (1|Species), family = gaussian(link = "identity"), df_sociality_matrix[df_sociality_matrix$k != 0,]))

### PGLS analysis with median k

df_sociality_matrix_median = 
  left_join(df_terminal_br_medians, sociality_matrix)

df_sociality_matrix_median$Sociality_score = as.numeric(df_sociality_matrix_median$Sociality_score)

cor.test(df_sociality_matrix_median$Sociality_score, df_sociality_matrix_median$k_median)

# t = 5.2093, df = 27, p-value = 1.736e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.4611630 0.8531318
# sample estimates:
#   cor 
# 0.7080028 

plot(df_sociality_matrix_median$Sociality_score, log2(df_sociality_matrix_median$k_median))

ggplot(df_sociality_matrix_median, aes(Sociality_score, k_median)) + 
  geom_point()

tree = read.nexus('/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/glmm/Ales_tree_modnames.txt')
tree$tip.label

df_sociality_matrix_median$Species = df_sociality_matrix_median$Sp_code

name.check(df_sociality_matrix_median, tree)


pglsModel = gls(k_median ~ Sociality_score, correlation = corPagel(1, phy = tree, fixed = F), df_sociality_matrix_median, method = 'ML')
# doesn't converge

summary(gls(k_median ~ Sociality_score, correlation = corBrownian(phy = tree), df_sociality_matrix_median, method = 'ML'))

summary(lm(k_median ~ Sociality_score, df_sociality_matrix_median))
