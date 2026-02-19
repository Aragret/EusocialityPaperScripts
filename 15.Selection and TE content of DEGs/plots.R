# setwd('/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/degs/')
setwd('M:/Share/PhD/eusociality_paper/analyses/degs/')

library(dplyr)
library(ggplot2)
library(emmeans)
library(cowplot)
library(lme4)

# import table with dNdS values
dnds_table = read.csv('dnds_table.csv') %>%
  mutate(gene_name = sub('-R.*', '', Node))

# import table with TE overlaps in 10kb upstream regions
TE_overlaps = read.csv('TE_10kbUp_overlaps_loI.csv',
                       header = FALSE) %>%
  mutate(gene_name = sub('-R.*', '', V1)) %>%
  select(gene_name, V2)

names(TE_overlaps) = c("gene_name", 'TE_overlaps_10kbUP')

# import table with HOG IDs
gene_HOGs = read.csv('S2_gene_to_hog.csv', sep = '\t', header = TRUE)
names(gene_HOGs) = c('gene_name', 'HOG')

# import table with RELAX results
relax = read.table('relax_table.tsv',
                   sep = '\t', header = T)
names(relax) = c('OG', 'LRT', 'p.value', 'k')

# import table with aBSREL results
absrel = read.table('absrel_significant_bifur',
                    header = F, col.names = 'OG') %>%
  mutate(Positive = 'Positive')

# import table with HOG and single-copy OG IDs
scOGs = read.table('SingleCopyHOGs',
                   header = F, sep = '\t') %>%
  mutate(HOG = sub('N0.', '', V1)) %>%
  select(HOG, V2)
names(scOGs) = c('HOG', 'OG')

# import table with DEGS and join with the previous tables
# DEG_pairwise&specificity_0.05.csv
deg_table = read.table('DEG_afterBC_per_caste_only_LFC_n005_FINAL.csv', header = TRUE, sep = ';') %>%
  left_join(dnds_table, by = 'gene_name') %>%
  left_join(TE_overlaps) %>% 
  left_join(gene_HOGs) %>%
  left_join(scOGs) %>%
  left_join(relax) %>%
  left_join(absrel) %>%
  mutate(BifurDev = as.factor(case_when(
    Species.x %in% c('Ncas', 'Kfla', 'Hsjo', 'Znev', 'PRsim', 'Bori', 'Bger', 'Cmer', 'Cpun') ~ 0,
    Species.x %in% c('Mdar', 'Rfla', 'Nluj', 'Apac', 'Mnat', 'Cges') ~ 1
  ))) %>%
  distinct()

# Specificity
#######
# Specificity

table(deg_table$Specificity)

table(deg_table[deg_table$Specificity == 'Wo',]$Species.x)
# Apac  Cges  Hsjo  Kfla  Mdar  Mnat  Ncas PRsim  Rfla 
# 605    24    16     4    79    68     2   238    31 

pdf('DEG_specificity_dnds.pdf')

ggplot(deg_table, aes(Specificity, w_ratios_Value)) +
  geom_jitter(aes(color = 'lightred', alpha = 0.7)) +
  stat_summary(fun.y=median, geom="point", shape=18, size=3) + theme(legend.position = "none")

ggplot(deg_table[deg_table$Specificity == 'Wo',], aes(Species.x, w_ratios_Value)) +
  geom_jitter(aes(color = 'lightred', alpha = 0.7)) +
  stat_summary(fun.y=median, geom="point", shape=18, size=3) + theme(legend.position = "none") +
  labs(title = 'dNdS ratios for worker-specific DEGs')

ggplot(deg_table[deg_table$Specificity == 'Wo',], aes(BifurDev, w_ratios_Value)) +
  geom_jitter(aes(color = 'lightred', alpha = 0.7)) +
  stat_summary(fun.y=median, geom="point", shape=18, size=3) + theme(legend.position = "none") +
  labs(title = 'dNdS ratios for worker-specific DEGs') +
  ylim(0, 2)

wilcox.test(deg_table[deg_table$Specificity == 'Wo' & deg_table$BifurDev == 0,]$w_ratios_Value,
            deg_table[deg_table$Specificity == 'Wo' & deg_table$BifurDev == 1,]$w_ratios_Value)
# p-value = 2.794e-09

# Compute the analysis of variance
res.aov <- aov(w_ratios_Value ~ Specificity, data = deg_table)
# Summary of the analysis
summary(res.aov)
# Specificity     9   57.4   6.372   48.67 <2e-16 ***

### TE overlap

ggplot(deg_table, aes(Specificity, TE_overlaps_10kbUP)) +
  geom_jitter(aes(color = 'lightred', alpha = 0.7)) +
  stat_summary(fun.y=median, geom="point", shape=18, size=3) + theme(legend.position = "none")

ggplot(deg_table[deg_table$Specificity == 'Wo',], aes(BifurDev, TE_overlaps_10kbUP)) +
  geom_jitter(aes(color = 'lightred', alpha = 0.7)) +
  stat_summary(fun.y=median, geom="point", shape=18, size=3) + theme(legend.position = "none")

wilcox.test(deg_table[deg_table$Specificity == 'Wo' & deg_table$BifurDev == 0,]$TE_overlaps_10kbUP,
            deg_table[deg_table$Specificity == 'Wo' & deg_table$BifurDev == 1,]$TE_overlaps_10kbUP)
# p-value < 2.2e-16

### RELAX

relax_deg = deg_table[!is.na(deg_table$k),]

ggplot(deg_table, aes(Specificity, k)) +
  geom_jitter(aes(color = 'lightred', alpha = 0.7)) +
  stat_summary(fun.y=median, geom="point", shape=18, size=3) + theme(legend.position = "none")

ggplot(deg_table, aes(Specificity, k)) +
  geom_jitter(aes(color = 'lightred', alpha = 0.7)) +
  stat_summary(fun.y=median, geom="point", shape=18, size=3) + theme(legend.position = "none") +
  ylim(0, 5)

ggplot(deg_table[deg_table$Specificity == 'Wo',], aes(BifurDev, k)) +
  geom_jitter(aes(color = 'lightred', alpha = 0.7)) +
  stat_summary(fun.y=median, geom="point", shape=18, size=3) + theme(legend.position = "none")

ggplot(deg_table[deg_table$Specificity == 'Wo',], aes(BifurDev, k)) +
  geom_jitter(aes(color = 'lightred', alpha = 0.7)) +
  stat_summary(fun.y=median, geom="point", shape=18, size=3) + theme(legend.position = "none") +
  ylim(0, 2)

wilcox.test(deg_table[deg_table$Specificity == 'Wo' & deg_table$BifurDev == 0,]$k,
            deg_table[deg_table$Specificity == 'Wo' & deg_table$BifurDev == 1,]$k)
# p-value = 0.236

dev.off()



########
### Specificity from 'fair' comparisons

deg_table_comp3 = deg_table %>%
  filter(Species.x %in% c('Mnat', 'Cges', 'Hsjo', 'Mdar', 'PRsim'))

# apac - no soldiers; kfla, ncas, rfla - many comparisons

ggplot(deg_table_comp3[deg_table_comp3$Specificity == 'Wo',], aes(BifurDev, w_ratios_Value)) +
  geom_jitter(aes(color = 'lightred', alpha = 0.7)) +
  stat_summary(fun.y=median, geom="point", shape=18, size=3) + theme(legend.position = "none") +
  labs(title = 'dNdS ratios for worker-specific DEGs')

wilcox.test(deg_table_comp3[deg_table_comp3$Specificity == 'Wo' & deg_table_comp3$BifurDev == 0,]$w_ratios_Value,
            deg_table_comp3[deg_table_comp3$Specificity == 'Wo' & deg_table_comp3$BifurDev == 1,]$w_ratios_Value)

# p-value = 0.01663

### aBSREL

absrel_deg = deg_table[!is.na(deg_table$Positive),]

table(absrel_deg$Specificity)

# Ad Al Er Ju Nm Ny Pr So Wo 
# 26 39 29 22 18 15 15  3 11 

table(absrel_deg$Specificity) / table(deg_table$Specificity) * 100

# Ad        Al        Er        Ju        Nm        Ny        Pr        So        Wo 
# 0.6535948 0.7500000 1.2144054 0.5062126 1.5503876 0.8108108 0.4509922 0.4109589 0.4552980 

### 3 branches

absrel_3branches = absrel_deg[absrel_deg$OG == 'OG0005794',] %>% distinct()

absrel_2branches = absrel_deg[absrel_deg$OG %in% c('OG0006166', 'OG0006468', 'OG0006558'),] %>% distinct()


### Repr vs workers


############## 

### Reproductive vs workers

ReprWo_deg = deg_table[, c('Species.x', 'gene_name', 'Al_Wo', 'Ne_Wo', 'Er_Wo', 'Nm_Wo', 'Pr_Wo', 'w_ratios_Value', 'TE_overlaps_10kbUP', 'k', 'BifurDev')] %>%
  filter(!(Species.x %in% c('Bori', 'BGER', 'Cmer', 'Cpun', 'Znev', 'Bger'))) %>%
  mutate(ReprWo = as.factor(case_when(
    Species.x == 'Mnat' & Al_Wo == 'Al' ~ 'Repr',
    Species.x == 'Mnat' & Al_Wo == 'Wo' ~ 'Wo',
    Species.x %in% c('Apac', 'Cges', 'PRsim') & Pr_Wo == 'Pr' ~ 'Repr',
    Species.x %in% c('Apac', 'Cges', 'PRsim') & Pr_Wo == 'Wo' ~ 'Wo',
    Species.x == 'Hsjo' & Ne_Wo == 'Ne' ~ 'Repr',
    Species.x == 'Hsjo' & Ne_Wo == 'Wo' ~ 'Wo',
    Species.x %in% c('Kfla', 'Ncas') & (Al_Wo == 'Al' | Ne_Wo == 'Ne') & Al_Wo != 'Wo' & Ne_Wo != 'Wo' ~ 'Repr',
    Species.x %in% c('Kfla', 'Ncas') & (Al_Wo == 'Wo' | Ne_Wo == 'Wo') & Al_Wo != 'Al' & Ne_Wo != 'Ne' ~ 'Wo',
    Species.x == 'Mdar' & Er_Wo == 'Er'  ~ 'Repr',
    Species.x == 'Mdar' & Er_Wo == 'Wo' ~ 'Wo',
    Species.x == 'Rfla' & Nm_Wo == 'Nm'  ~ 'Repr',
    Species.x == 'Rfla' & Nm_Wo == 'Wo' ~ 'Wo',
    .default = 'None'
  ))) %>%
  distinct()

ReprWo_deg$Species.x = factor(ReprWo_deg$Species.x, levels = c('Mdar', 'Hsjo', 'Kfla', 'Ncas', 'PRsim', 'Rfla', 'Cges', 'Mnat', 'Apac'))

write.csv(ReprWo_deg, 'ReprWo_deg.csv', quote = F, row.names = F)

table(ReprWo_deg[ReprWo_deg$ReprWo == 'Wo',]$Species.x)
# Apac  Cges  Hsjo  Kfla  Mdar  Mnat  Ncas PRsim  Rfla 
# 603   671    52  1870  1471  1717  2936  1066  1344 

table(ReprWo_deg[ReprWo_deg$ReprWo == 'Repr',]$Species.x)
# Apac  Cges  Hsjo  Kfla  Mdar  Mnat  Ncas PRsim  Rfla 
# 686   875    88  2151  1493  1708  2895  1470  1431 

# worker-biased and -specific in Apac
wilcox.test(ReprWo_deg[ReprWo_deg$ReprWo == 'Wo' & ReprWo_deg$Species.x == 'Apac',]$w_ratios_Value,
            deg_table[deg_table$Specificity == 'Wo' & ReprWo_deg$Species.x == 'Apac',]$w_ratios_Value)
# p-value = 0.4338

ggplot(ReprWo_deg[ReprWo_deg$ReprWo == 'Wo',], aes(Species.x, w_ratios_Value)) +
  geom_jitter(aes(color = 'lightred', alpha = 0.7)) +
  stat_summary(fun.y=median, geom="point", shape=18, size=3) + theme(legend.position = "none") +
  labs(title = 'dNdS ratios for worker-biased DEGs')


ggplot(ReprWo_deg[ReprWo_deg$ReprWo == 'Wo',], aes(BifurDev, w_ratios_Value)) +
  geom_boxplot() +
  ylim(0, 1)

ggplot(ReprWo_deg, aes(Species.x, w_ratios_Value, fill = ReprWo)) +
  geom_boxplot()
  # geom_jitter(aes(color = BifurDev, alpha = 0.7)) +
  # stat_summary(fun.y=median, geom="point", shape=18, size=3)

ReprWo_w_species = ggplot(ReprWo_deg, aes(Species.x, w_ratios_Value, fill = ReprWo)) +
  geom_boxplot() +
  ylim(0, 1) 

ReprWo_w_species

ggplot(ReprWo_deg, aes(ReprWo, w_ratios_Value, fill = BifurDev)) +
  geom_boxplot() +
  ylim(0, 1)

ggplot(ReprWo_deg, aes(Species.x, TE_overlaps_10kbUP, fill = ReprWo)) +
  geom_boxplot() +
  ylim(0, 20)

ggplot(ReprWo_deg, aes(ReprWo, TE_overlaps_10kbUP, fill = BifurDev)) +
  geom_boxplot() +
  ylim(0, 20)

ggplot(ReprWo_deg, aes(Species.x, k, fill = ReprWo)) +
  geom_boxplot() +
  ylim(0, 5)

ggplot(ReprWo_deg, aes(ReprWo, k, fill = BifurDev)) +
  geom_boxplot() +
  ylim(0, 2)


### without cockroaches
######


# CA
mod= glm(w_ratios_Value ~ BifurDev * ReprWo, data = ReprWo_deg)
anova(mod)
summary(mod)
str(ReprWo_deg)
em = emmeans(mod, pairwise ~ BifurDev + ReprWo)
pairs(em , simple = "each")
# end

# $`comparisons of simple contrasts for BifurDev`
# ReprWo = None:
#   contrast              estimate      SE    df t.ratio p.value
# BifurDev0 - BifurDev1 -0.03257 0.00312 94447 -10.455  <.0001
# 
# ReprWo = Repr:
#   contrast              estimate      SE    df t.ratio p.value
# BifurDev0 - BifurDev1 -0.09325 0.00823 94447 -11.329  <.0001
# 
# ReprWo = Wo:
#   contrast              estimate      SE    df t.ratio p.value
# BifurDev0 - BifurDev1  0.00126 0.00879 94447   0.144  0.8858
# 
# 
# $`comparisons of simple contrasts for ReprWo`
# BifurDev = 0:
#   contrast    estimate      SE    df t.ratio p.value
# None - Repr   0.0834 0.00620 94447  13.450  <.0001
# None - Wo     0.0289 0.00668 94447   4.320  <.0001
# Repr - Wo    -0.0546 0.00849 94447  -6.424  <.0001
# 
# BifurDev = 1:
#   contrast    estimate      SE    df t.ratio p.value
# None - Repr   0.0228 0.00624 94447   3.646  0.0008
# None - Wo     0.0627 0.00650 94447   9.643  <.0001
# Repr - Wo     0.0399 0.00853 94447   4.681  <.0001
# 
# P value adjustment: tukey method for comparing a family of 3 estimates 

modTE = glm(TE_overlaps_10kbUP ~ BifurDev * ReprWo, data = ReprWo_deg)
summary(modTE)

emTE = emmeans(modTE, pairwise ~ BifurDev + ReprWo)
pairs(emTE, simple = 'each')

# ReprWo = None:
#   contrast              estimate     SE     df t.ratio p.value
# BifurDev0 - BifurDev1   -0.338 0.0271 149085 -12.472  <.0001
# 
# ReprWo = Repr:
#   contrast              estimate     SE     df t.ratio p.value
# BifurDev0 - BifurDev1    1.178 0.0838 149085  14.064  <.0001
# 
# ReprWo = Wo:
#   contrast              estimate     SE     df t.ratio p.value
# BifurDev0 - BifurDev1    0.552 0.0875 149085   6.316  <.0001
# 
# 
# $`comparisons of simple contrasts for ReprWo`
# BifurDev = 0:
#   contrast    estimate     SE     df t.ratio p.value
# None - Repr   -1.023 0.0618 149085 -16.572  <.0001
# None - Wo     -0.510 0.0648 149085  -7.870  <.0001
# Repr - Wo      0.513 0.0847 149085   6.053  <.0001
# 
# BifurDev = 1:
#   contrast    estimate     SE     df t.ratio p.value
# None - Repr    0.493 0.0627 149085   7.852  <.0001
# None - Wo      0.380 0.0646 149085   5.879  <.0001
# Repr - Wo     -0.113 0.0865 149085  -1.302  0.3938
# 
# P value adjustment: tukey method for comparing a family of 3 estimates 

modK = glm(k ~ BifurDev * ReprWo, data = ReprWo_deg)
summary(modK)

emK = emmeans(modK, pairwise ~ BifurDev + ReprWo)
pairs(emK, simple = 'each')

# $`comparisons of simple contrasts for BifurDev`
# ReprWo = None:
#   contrast              estimate    SE   df t.ratio p.value
# BifurDev0 - BifurDev1   -0.124 0.150 7293  -0.824  0.4101
# 
# ReprWo = Repr:
#   contrast              estimate    SE   df t.ratio p.value
# BifurDev0 - BifurDev1    0.147 0.378 7293   0.389  0.6972
# 
# ReprWo = Wo:
#   contrast              estimate    SE   df t.ratio p.value
# BifurDev0 - BifurDev1    0.707 0.419 7293   1.687  0.0916
# 
# 
# $`comparisons of simple contrasts for ReprWo`
# BifurDev = 0:
#   contrast    estimate    SE   df t.ratio p.value
# None - Repr   0.0786 0.298 7293   0.264  0.9624
# None - Wo    -0.6429 0.305 7293  -2.111  0.0877
# Repr - Wo    -0.7214 0.395 7293  -1.827  0.1608
# 
# BifurDev = 1:
#   contrast    estimate    SE   df t.ratio p.value
# None - Repr   0.3492 0.276 7293   1.264  0.4156
# None - Wo     0.1877 0.325 7293   0.578  0.8317
# Repr - Wo    -0.1614 0.403 7293  -0.401  0.9154
# 
# P value adjustment: tukey method for comparing a family of 3 estimates 


### MANOVA
##############################################################################

summary(aov(w_ratios_Value ~ BifurDev + ReprWo, data = ReprWo_deg))
# Df Sum Sq Mean Sq F value Pr(>F)    
# BifurDev        1   12.0  11.971 104.224 <2e-16 ***
#   ReprWo          1    0.3   0.279   2.428  0.119    
# Residuals   19755 2269.0   0.115                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 129394 Beobachtungen als fehlend gelöscht

summary(aov(w_ratios_Value ~ BifurDev * ReprWo, data = ReprWo_deg))
# Df Sum Sq Mean Sq F value Pr(>F)    
# BifurDev            1   12.0  11.971 104.725 <2e-16 ***
#   ReprWo              1    0.3   0.279   2.439  0.118    
# BifurDev:ReprWo     1   11.0  10.979  96.047 <2e-16 ***
#   Residuals       19754 2258.0   0.114                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 129394 Beobachtungen als fehlend gelöscht

summary(aov(TE_overlaps_10kbUP ~ BifurDev * ReprWo, data = ReprWo_deg))
# Df Sum Sq Mean Sq F value   Pr(>F)    
# BifurDev            1   4756    4756  241.95  < 2e-16 ***
#   ReprWo              1    261     261   13.28 0.000269 ***
#   BifurDev:ReprWo     1    598     598   30.44 3.47e-08 ***
#   Residuals       24518 481935      20                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 124630 Beobachtungen als fehlend gelöscht

summary(aov(k ~ BifurDev * ReprWo, data = ReprWo_deg))
# Df Sum Sq Mean Sq F value Pr(>F)
# BifurDev           1     75   75.05   2.362  0.125
# ReprWo             1     79   78.91   2.484  0.115
# BifurDev:ReprWo    1     31   30.93   0.974  0.324
# Residuals       1601  50867   31.77               
# 147547 Beobachtungen als fehlend gelöscht

####

pdf('ReprWo.pdf')

ggplot(ReprWo_deg, aes(BifurDev, w_ratios_Value, fill = ReprWo)) +
  geom_boxplot() + 
  ylim(0, 1)

ggplot(ReprWo_deg, aes(BifurDev, TE_overlaps_10kbUP, fill = ReprWo)) +
  geom_boxplot() + 
  ylim(0, 30)

ggplot(ReprWo_deg, aes(BifurDev, k, fill = ReprWo)) +
  geom_boxplot() + 
  ylim(0, 1.5)

dev.off()






###############################################################
### Cockroaches

deg_table_roaches = deg_table[deg_table$Species.x %in% 
                                c('Bori', 'Bger', 'Cmer', 'Cpun'),]

deg_table_roaches[is.na(deg_table_roaches$Specificity),]$Specificity = 'None'

ggplot(deg_table_roaches, aes(Specificity, TE_overlaps_10kbUP)) +
  geom_boxplot() +
  ylim(0, 25)

roachTE = ggplot(deg_table_roaches[deg_table_roaches$Species.x %in% c('Bori', 'Cmer'),], aes(Species.x, TE_overlaps_10kbUP, fill = Specificity)) +
  geom_boxplot(width=0.15*length(unique(deg_table_roaches[deg_table_roaches$Species.x %in% c('Bori', 'Cmer'),]$Species.x))) +
  ylim(0, 25) + theme_classic()

### combine with termites

ReprWo_deg$ReprWo = factor(ReprWo_deg$ReprWo, levels = c('Repr', 'Wo', 'None'))

linTE = ggplot(ReprWo_deg[ReprWo_deg$BifurDev == 0,], aes(Species.x, TE_overlaps_10kbUP, fill = ReprWo)) +
  geom_boxplot(width=0.1*length(unique(ReprWo_deg[ReprWo_deg$BifurDev == 0,]$Species.x))) +
  ylim(0, 25) + 
  theme_classic() + guides(fill = 'none')

bifurTE = ggplot(ReprWo_deg[ReprWo_deg$BifurDev == 1,], aes(Species.x, TE_overlaps_10kbUP, fill = ReprWo)) +
  geom_boxplot(width=0.15*length(unique(ReprWo_deg[ReprWo_deg$BifurDev == 1,]$Species.x))) +
  ylim(0, 25) + theme_classic()

TE_plots = plot_grid(roachTE, linTE, bifurTE, nrow = 1)

save_plot('TE_DEGs.pdf', TE_plots, base_width = 11)
# 
# modTEroaches = glm(TE_overlaps_10kbUP ~ Specificity, data = deg_table_roaches)
# summary(modTEroaches)
# 
# emTE = emmeans(modTEroaches, pairwise ~ Specificity)
# pairs(emTE, simple = 'each')

# contrast estimate   SE   df t.ratio p.value
# Ad - Ju   -0.0275 0.21 1680  -0.131  0.8957

ggplot(deg_table_roaches, aes(Specificity, w_ratios_Value)) +
  geom_boxplot() +
  ylim(0, 1.5)

roachesW = ggplot(deg_table_roaches, aes(Species.x, w_ratios_Value, fill = Specificity)) +
  geom_boxplot(width=0.15*length(unique(deg_table_roaches$Species.x))) +
  ylim(0, 1) + theme_classic()

linW = ggplot(ReprWo_deg[ReprWo_deg$BifurDev == 0,], aes(Species.x, w_ratios_Value, fill = ReprWo)) +
  geom_boxplot(width=0.1*length(unique(ReprWo_deg[ReprWo_deg$BifurDev == 0,]$Species.x))) +
  ylim(0, 1) + 
  theme_classic() + guides(fill = 'none')

bifurW = ggplot(ReprWo_deg[ReprWo_deg$BifurDev == 1,], aes(Species.x, w_ratios_Value, fill = ReprWo)) +
  geom_boxplot(width=0.15*length(unique(ReprWo_deg[ReprWo_deg$BifurDev == 1,]$Species.x))) +
  ylim(0, 1) + theme_classic()

W_plots = plot_grid(roachesW, linW, bifurW, nrow = 1)

save_plot('dNdS_DEGs.pdf', W_plots, base_width = 11)

# modroaches = glm(w_ratios_Value ~ Specificity, data = deg_table_roaches)
# summary(modroaches)
# 
# emmod = emmeans(modroaches, pairwise ~ Specificity)
# pairs(emmod, simple = 'each')

# contrast estimate     SE   df t.ratio p.value
# Ad - Ju     0.057 0.0119 3211   4.808  <.0001

ggplot(deg_table_roaches, aes(Specificity, k)) +
  geom_boxplot() +
  ylim(0, 2)

modKroaches = glm(k ~ Specificity, data = deg_table_roaches)
summary(modKroaches)

emk = emmeans(modKroaches, pairwise ~ Specificity)
pairs(emk, simple = 'each')

# contrast estimate  SE  df t.ratio p.value
# Ad - Ju      -1.8 0.9 259  -2.000  0.0466

########################################
### combine termite and roaches stats

deg_roaches = deg_table_roaches %>%
  mutate(ReprWo = case_when(
    Specificity == 'None' ~ 'None',
    Specificity == 'Ju' ~ 'Wo',
    Specificity == 'Ad' ~ 'Repr'
  ),
  # BifurDev = '2') %>%
  BifurDev = case_when(
    Species.x %in% c('Bger', 'Bori') ~ '2',
    Species.x %in% c('Cmer', 'Cpun') ~ '3'
  )) %>%
  select(Species.x, gene_name, ReprWo, w_ratios_Value, TE_overlaps_10kbUP,
         k, BifurDev)

# deg_roaches$Termite = 0

ReprWo_all = ReprWo_deg %>%
  select(Species.x, gene_name, ReprWo, w_ratios_Value, TE_overlaps_10kbUP,
         k, BifurDev) %>%
  mutate(BifurDev = as.character(BifurDev)) %>%
#   mutate(Termite = 1) %>%
  bind_rows(deg_roaches)


ReprWo_all$Species.x = as.factor(ReprWo_all$Species.x)

modTE = glmer(TE_overlaps_10kbUP ~ BifurDev * ReprWo + (1|Species.x), data = ReprWo_all)
summary(modTE)

emTE = emmeans(modTE, pairwise ~ ReprWo + BifurDev)
pairs(emTE, simple = 'each')

# BifurDev = 0:
#   contrast    estimate     SE  df z.ratio p.value
# None - Repr  -0.4952 0.0587 Inf  -8.433  <.0001
# None - Wo     0.1457 0.0618 Inf   2.356  0.0484
# Repr - Wo     0.6409 0.0782 Inf   8.200  <.0001
# 
# BifurDev = 1:
#   contrast    estimate     SE  df z.ratio p.value
# None - Repr   0.0803 0.0582 Inf   1.380  0.3517
# None - Wo    -0.0321 0.0601 Inf  -0.534  0.8546
# Repr - Wo    -0.1124 0.0797 Inf  -1.410  0.3358
# 
# BifurDev = 2:
#   contrast    estimate     SE  df z.ratio p.value
# None - Repr  -0.6304 0.1820 Inf  -3.457  0.0016
# None - Wo    -0.4827 0.1330 Inf  -3.617  0.0009
# Repr - Wo     0.1477 0.2240 Inf   0.660  0.7865
# 
# Degrees-of-freedom method: asymptotic 
# P value adjustment: tukey method for comparing a family of 3 estimates 
# 
# $`comparisons of simple contrasts for BifurDev`
# ReprWo = None:
#   contrast              estimate   SE  df z.ratio p.value
# BifurDev0 - BifurDev1   -0.126 1.20 Inf  -0.105  0.9940
# BifurDev0 - BifurDev2   -2.179 1.55 Inf  -1.404  0.3387
# BifurDev1 - BifurDev2   -2.053 1.50 Inf  -1.370  0.3569
# 
# ReprWo = Repr:
#   contrast              estimate   SE  df z.ratio p.value
# BifurDev0 - BifurDev1    0.450 1.20 Inf   0.374  0.9260
# BifurDev0 - BifurDev2   -2.314 1.56 Inf  -1.481  0.3002
# BifurDev1 - BifurDev2   -2.764 1.51 Inf  -1.829  0.1600
# 
# ReprWo = Wo:
#   contrast              estimate   SE  df z.ratio p.value
# BifurDev0 - BifurDev1   -0.304 1.20 Inf  -0.252  0.9656
# BifurDev0 - BifurDev2   -2.807 1.56 Inf  -1.802  0.1690
# BifurDev1 - BifurDev2   -2.504 1.51 Inf  -1.663  0.2196
# 
# Degrees-of-freedom method: asymptotic 
# P value adjustment: tukey method for comparing a family of 3 estimates 



ReprWo_all %>%
  group_by(BifurDev) %>%
  mutate(Mean = mean(TE_overlaps_10kbUP))

by(ReprWo_all, ReprWo_all$BifurDev, summary)



modW = glmer(w_ratios_Value ~ BifurDev * ReprWo + (1|Species.x), data = ReprWo_all)
summary(modW)

emW = emmeans(modW, pairwise ~ BifurDev + ReprWo)
pairs(emW, simple = 'each')

# $`comparisons of simple contrasts for BifurDev`
# ReprWo = None:
#   contrast              estimate     SE  df z.ratio p.value
# BifurDev0 - BifurDev1 -0.03084 0.0353 Inf  -0.873  0.6572
# BifurDev0 - BifurDev2  0.05192 0.0372 Inf   1.395  0.3436
# BifurDev1 - BifurDev2  0.08277 0.0353 Inf   2.344  0.0500
# 
# ReprWo = Repr:
#   contrast              estimate     SE  df z.ratio p.value
# BifurDev0 - BifurDev1 -0.09799 0.0362 Inf  -2.708  0.0186
# BifurDev0 - BifurDev2 -0.01344 0.0393 Inf  -0.342  0.9375
# BifurDev1 - BifurDev2  0.08455 0.0375 Inf   2.254  0.0624
# 
# ReprWo = Wo:
#   contrast              estimate     SE  df z.ratio p.value
# BifurDev0 - BifurDev1 -0.00513 0.0363 Inf  -0.141  0.9891
# BifurDev0 - BifurDev2  0.08481 0.0390 Inf   2.174  0.0756
# BifurDev1 - BifurDev2  0.08994 0.0372 Inf   2.420  0.0411
# 
# Degrees-of-freedom method: asymptotic 
# P value adjustment: tukey method for comparing a family of 3 estimates 
# 
# $`comparisons of simple contrasts for ReprWo`
# BifurDev = 0:
#   contrast    estimate      SE  df z.ratio p.value
# None - Repr   0.0787 0.00656 Inf  11.992  <.0001
# None - Wo     0.0252 0.00706 Inf   3.573  0.0010
# Repr - Wo    -0.0535 0.00865 Inf  -6.183  <.0001
# 
# BifurDev = 1:
#   contrast    estimate      SE  df z.ratio p.value
# None - Repr   0.0116 0.00640 Inf   1.808  0.1670
# None - Wo     0.0509 0.00668 Inf   7.628  <.0001
# Repr - Wo     0.0394 0.00868 Inf   4.538  <.0001
# 
# BifurDev = 2:
#   contrast    estimate      SE  df z.ratio p.value
# None - Repr   0.0133 0.01180 Inf   1.133  0.4935
# None - Wo     0.0581 0.01050 Inf   5.552  <.0001
# Repr - Wo     0.0448 0.01530 Inf   2.920  0.0098
# 
# Degrees-of-freedom method: asymptotic 
# P value adjustment: tukey method for comparing a family of 3 estimates 



by(ReprWo_all, ReprWo_all[, c('BifurDev', 'ReprWo')], summary)
