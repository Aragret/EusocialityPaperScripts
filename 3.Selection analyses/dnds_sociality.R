library(dplyr)
library(ape)
library(lme4)
library(geiger)
library(nlme)
library(phytools)


setwd('/global/scratch2/amikhail/termite_genomes/eusociality_paper/codeml_allHOGs_results/')

dnds_file = read.csv('dNdSw_table.csv')
dnds_file_2 = read.csv('dNdSw_table_2.csv')

dnds_table = rbind(dnds_file, dnds_file_2) %>%
  filter(!grepl('InternalNode', Node)) %>%
  mutate(Species = substr(Node, 1, 4)) %>%
  filter(w_ratios_Value < 10) %>%
  filter(ds_Value < 3)

tree = read.nexus('/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/glmm/Ales_tree_modnames.txt')

sociality_matrix = read.csv('Copy_of_Matrix_of_ecological traits_13.06.2024 - Formated_matrix.csv', dec = '.') %>%
  select(Ontogeny, Nesting_type, Colony_size_log10_lower_range, Sociality_score, Sp_code) %>%
  filter(Sp_code != '')

sociality_matrix$Species = sociality_matrix$Sp_code

df = merge(dnds_table, sociality_matrix, by.x = 'Species', by.y = 'Sp_code')

df_medians = df %>%
  group_by(Species) %>%
  summarize(w_ratios_median = median(w_ratios_Value)) %>%
  left_join(sociality_matrix)

summary(glmer(w_ratios_Value ~ Sociality_score + (1|Orthogroup) + (1|Species), family=Gamma(link = "log"), df))

# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: Gamma  ( log )
# Formula: w_ratios_Value ~ Sociality_score + (1 | Orthogroup) + (1 | Species)
# Data: df
# 
# AIC       BIC    logLik  deviance  df.resid 
# -317600.3 -317547.3  158805.2 -317610.3    300964 
# 
# Scaled residuals: 
#   Min     1Q Median     3Q    Max 
# -0.849 -0.470 -0.175  0.185 35.322 
# 
# Random effects:
#   Groups     Name        Variance Std.Dev.
# Orthogroup (Intercept) 1.27918  1.1310  
# Species    (Intercept) 0.05577  0.2362  
# Residual               1.38879  1.1785  
# Number of obs: 300969, groups:  Orthogroup, 13522; Species, 29
# 
# Fixed effects:
#   Estimate Std. Error t value Pr(>|z|)    
# (Intercept)     -1.48305    0.05178 -28.641  < 2e-16 ***
#   Sociality_score -0.09622    0.03057  -3.148  0.00164 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr)
# Socilty_scr -0.676

# summary(glmer(w_ratios_Value ~ Colony_size_log10_lower_range + (1|Orthogroup) + (1|Species), family=Gamma(link = "log"), df))


name.check(df_medians, tree)

plot(df_medians$Sociality_score, df_medians$w_ratios_median)

pglsModel = gls(Sociality_score ~ w_ratios_median, correlation = corPagel(1, phy = tree, fixed = F), df_medians, method = 'ML')
summary(pglsModel)

# Model: Sociality_score ~ w_ratios_median 
# Data: df_medians 
# AIC      BIC    logLik
# 81.48024 86.94942 -36.74012
# 
# Correlation Structure: corPagel
# Formula: ~1 
# Parameter estimate(s):
#   lambda 
# 1.051206 
# 
# Coefficients:
#   Value Std.Error   t-value p-value
# (Intercept)       3.434893  1.021790  3.361643  0.0023
# w_ratios_median -10.470120  2.460818 -4.254732  0.0002
# 
# Correlation: 
#   (Intr)
# w_ratios_median -0.45 
# 
# Standardized residuals:
#   Min         Q1        Med         Q3        Max 
# -0.6400247 -0.4856471 -0.3695483 -0.1801084  2.1350315 
# 
# Residual standard error: 1.763368 
# Degrees of freedom: 29 total; 27 residual

summary(gls(w_ratios_median ~ Sociality_score, correlation = corBrownian(phy = tree), df_medians, method = 'ML'))

# corPagel doesn't work

# Generalized least squares fit by maximum likelihood
# Model: w_ratios_median ~ Sociality_score 
# Data: df_medians 
# AIC       BIC   logLik
# -90.04622 -85.94433 48.02311
# 
# Correlation Structure: corBrownian
# Formula: ~1 
# Parameter estimate(s):
#   numeric(0)
# 
# Coefficients:
#   Value  Std.Error   t-value p-value
# (Intercept)      0.23601868 0.04807782  4.909097  0.0000
# Sociality_score -0.02776954 0.00759713 -3.655269  0.0011
# 
# Correlation: 
#   (Intr)
# Sociality_score -0.227
# 
# Standardized residuals:
#   Min            Q1           Med            Q3           Max 
# -1.2642253082 -0.4586341026 -0.2347224093  0.0006332402  0.7699862537 
# 
# Residual standard error: 0.08800876 
# Degrees of freedom: 29 total; 27 residual

null_model = gls(Sociality_score ~ 1, correlation = corPagel(1, phy = tree, fixed = F), df_medians, method = 'ML')

coef(pglsModel)
coef(null_model)

anova(null_model, pglsModel)

