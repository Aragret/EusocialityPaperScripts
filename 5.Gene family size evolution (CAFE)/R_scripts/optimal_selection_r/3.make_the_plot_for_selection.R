################################################################################
#
#             Title: Optimal model selection 
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
Here we have already run the tree wit the optimal number of k and l. We check 
if the model converge or not. If not, we test a simpler one.

We first test k=2, l=4. The model do not converge

'


# ----- [0.2] Dependancies: ----------------------------------------------------

'The following libraries are required: ggplot2'
'The following input files are required:
optimal_selection.txt
The output files are:
opt_selection_plot_lambda_*
opt_selection_plot_AIC_*
AIC_table_KL_models.txt
Optimal_selection_withAIC_values.txt
'

# ------------------------------------------------------------------------------
# ------ [1] Introduction ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- [1.1] Load packages ----------------------------------------------------

library(ggplot2)

# ----- [1.2] Set directory ----------------------------------------------------



pa = "C:/Users/caumont/Documents/WORK/Dino_McMahon/Work/Lab work and projects/"
setwd(
  paste0(
    pa, "10. Eusociality/CAFE/4. Selection"
  )
)


# ----- [1.3] Load data --------------------------------------------------------

data <-read.table(file= 'optimal_selection.txt',
                  header=T,
                  dec=".",
                  stringsAsFactors = T,
                  fill = T, numerals = "no.loss")
str(data)
data$AIC = 2*(data$no_rate + data$no_cluster) - 2*(-data$lnL)
colnames(data) = c("model", "no_rate","rate", "no_cluster", "L", "li", "likelihood", "lambda", "AIC")
datasub = subset(data, data$L %in% c("L1","L2","L3","L4"))


# ------- Function plot and save -----------------------------------------------


plot_save = function(data, string){
  ggplot(data = data, aes(li,lambda))+
    geom_boxplot(aes(fill=L), outlier.shape = NA, show.legend = F )+
    facet_grid(rows = vars(rate), cols = vars(L))+
    geom_jitter(alpha=0.3, height = 0)+
    ylim(0,9)
  
  ggsave(filename = paste0("opt_selection_plot_lambda_",string,".jpg"),
         width = 14,height = 8,units = "in", 
         dpi = 600, limitsize = F)
  
  ggplot(data = data, aes(li,AIC))+
    stat_boxplot(aes(fill=L) ,outlier.shape = NA,show.legend = F )+
    facet_grid(rows = vars(rate), cols = vars(L))+
    geom_jitter(alpha=0.3, height = 0)+
    ylab("AIC")+
    ylim(600000,900000)
  
  ggsave(filename = paste0("opt_selection_plot_AIC_",string,".jpg"),
         width = 14,height = 8,units = "in", 
         dpi = 600, limitsize = F)
}


# ------------------------------------------------------------------------------

plot_save(data = data, "Eusociality")

plot_save(data = datasub, "Eusociality_subdataset")

# ---------------------- Median AIC --------------------------------------------

model_list = as.character(unique(data$model))
AIC_table = data.frame(list("Model" = model_list, "AIC" = "NA"))

for (model in model_list){
  Med = median(data$AIC[data$model==model])
  AIC_table$AIC[AIC_table$Model==model] = round(Med,2)
}
AIC_table$AIC = as.numeric(AIC_table$AIC)
AIC_table

write.table(AIC_table, file = "AIC_table_KL_models.txt",quote = F,
            dec = ",", sep = "\t",row.names = F)
write.table(data,file = "Optimal_selection_withAIC_values.txt",quote = F,
            dec = ",", sep = "\t",row.names = F)
