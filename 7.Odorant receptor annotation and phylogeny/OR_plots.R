# install.packages('itol.toolkit')
# install.packages('ape')
# install.packages('readxl')

library(itol.toolkit)
library(ape)
library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(tibble)
library(stringr)

tree = read.tree('M:/Share/ORs_rooted.txt')
# tree$tip.label = sub('-R.*', '', tree$tip.label)

dataset = as.data.frame(tree$tip.label)
dataset$Species = unlist(lapply(dataset$`tree$tip.label`, substr, 0, 4))

# unit_12 <- create_unit(data = dataset, 
#                        key = "E012_style_1",
#                        type = "DATASET_STYLE", 
#                        subtype = "branch",
#                        position = "clade",
#                        size_factor = 5,
#                        tree = tree)
# write_unit(unit_12, file = 'M:/Share/all_ORs_einsi_outgroup_trimal.msa.treefile_metadata')

dataset_groups = dataset %>%
  mutate(
    Group = case_when(
      Species %in% c('BGER', 'Bori', 'KAJ4', 'Dpun', 'Bger', 'Pame', 'GR_B') ~ 'Cockroaches',
      Species %in% c('Cmer', 'Cpun') ~ 'Cryptocercus',
      # Species == 'Mdar' ~ 'Mastotermes',
      # Species == 'Dlon' ~ 'Dolichorhinotermes',
      Species %in% c('Hsjo', 'Znev', 'Ncas', 'Kfla', 'PRsi', 'Cbre', 'Isch') ~ 'TermitesLinear',
      Species %in% c('Mdar', 'Dlon', 'Rfla', 'Hten', 'Ctes', 'Cges', 'Ofor', 'Mnat', 'Ssph', 'Aunk', 'Apac', 'Nluj', 'Csp4', 'Ntar', 'Pred', 'Munk') ~
        'TermitesBifurcated'
    )
  ) %>%
  select(-Species)



Group_color <- tribble(
  ~Group, ~Color,
  "Cockroaches", "#dededeff",
  "Cryptocercus", "#8a8a8aff",
  "TermitesLinear", "#000000ff",
  "TermitesBifurcated", "#12c1edff"
)

# Group_color <- tribble(
#   ~Group, ~Color,
#   "Cockroaches", "#edf8b1",
#   "Cryptocercus", "#7fcdbb",
#   "Termites", "#2c7fb8"
# )

Group_data = left_join(dataset_groups, Group_color)



unit_12_groups <- create_unit(data = Group_data, 
                             key = "color_manual_1",
                             type = "DATASET_STYLE", 
                             subtype = "branch",
                             position = "clade",
                             size_factor = 14,
                             tree = tree)
write_unit(unit_12_groups, file = 'M:/Share/all_ORs_einsi_outgroup_trimal.msa.treefile_metadata_groups')
