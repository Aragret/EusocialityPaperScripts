library(dplyr)
library(ggplot2)
library(cowplot)

# setwd("~/EusocialityPaper/relax/")
setwd('M:/Share/eusociality_paper/selection/relax/')

df = read.csv("bootstrap_relax.csv") %>%
  group_by(Bootstrap) %>%
  mutate(p.adj = p.adjust(p.value, method = 'BH'))

df_significant = df[df$p.adj < 0.05,]

df_sign_summary = df_significant %>%
  group_by(Bootstrap) %>%
  summarise()

df_sign_relax = df_significant[df_significant$k < 1,]
df_sign_intens = df_significant[df_significant$k > 1,]

relax = as.data.frame(table(df_sign_relax$Bootstrap))

relax$Bifur = 10

par(mfrow = c(1, 2))

first = ggplot(relax, aes("Relaxed", as.numeric(Freq))) +
  geom_boxplot() +
  geom_point(aes(y = Bifur), color="#12c1edff", size = 3) + 
  theme_minimal() + ylab('Frequency') + 
  coord_cartesian(ylim = c(0, 40))

intens = as.data.frame(table(df_sign_intens$Bootstrap))

intens$Bifur = 5

second = ggplot(intens, aes("Intensified", as.numeric(Freq))) +
  geom_boxplot() +
  geom_point(aes(y = Bifur), color="#12c1edff", size = 3) +
  theme_minimal() + ylab('') +
  coord_cartesian(ylim = c(0, 60))

distr_relax = ecdf(relax$Freq)
distr_relax(10)*100 # [1] 61.61616

distr_intens = ecdf(intens$Freq)
distr_intens(5)*100 # [1] 13.26531

save_plot("Bootstrap_RELAX.pdf", plot_grid(first, second), nrow = 1)
