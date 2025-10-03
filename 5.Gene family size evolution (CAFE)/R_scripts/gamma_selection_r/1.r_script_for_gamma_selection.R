################################################################################
#
#             Title: Gamma selection 
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
This script is used to plot the variation in lamda at different gamma rate to
check for convergence and then find the model with the smallest likelihood.
'


# ----- [0.2] Dependancies: ----------------------------------------------------

'The following libraries are required: ggplot2'
'The following input files are required: 
gamma_selection.txt'
'The output files are: 
gamma_selection_eusoc_lambda.tiff
gamma_selection_plot_AIC.tiff
gamma_selection_plot_med_lambda.tiff
gamma_selection_plot_med_lnL.tiff
'

# ------------------------------------------------------------------------------
# ------ [1] Introduction ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- [1.1] Load packages ----------------------------------------------------

library(ggplot2)

# ----- [1.2] Set directory ----------------------------------------------------

pa = "~/WORK/Dino_McMahon/Work/Lab work and projects/"
setwd(
  paste0(
    pa, "10. Eusociality/CAFE/4. Selection"
  )
)

# ----- [1.3] convergence function ---------------------------------------------

'Function which calculate the median value of lambda for each gamma rate
categories based on 50 runs. Lambda values for each run are selected into 
another dataset if they are close to the median'
keep_median_lambda_values = function(dataset){
  lnl_list = c()
  lambda_list = c()
  run_list = c()
  rate_list = c()
  for (M in levels(dataset$rate)){
    B_M1_values = dataset[dataset$rate %in% c(M),3]
    B_M1_rounds = round(B_M1_values,3)
    B_M1_rndmed = median(B_M1_rounds)
    B_M1_median = median(B_M1_values)
    B_M1_select = B_M1_values[B_M1_rounds==B_M1_rndmed]
    len_select = length(B_M1_select)
    rate_list = append(rate_list, rep(M, len_select))
    lambda_list = append(lambda_list, B_M1_select)
    lnl_list = append(lnl_list, dataset$likelihood[dataset$rate==M & dataset$lambda %in% B_M1_select])
  }
  new_dataset = data.frame(cbind(rate_list,lambda_list,lnl_list))
  colnames(new_dataset)= c("rate", "lambda", "likelihood")
  
  new_dataset$rate = as.factor(new_dataset$rate)
  new_dataset$lambda = as.numeric(new_dataset$lambda)
  new_dataset$likelihood = as.numeric(new_dataset$likelihood)
  
  
  new_dataset$rate=ordered(new_dataset$rate, levels=c("k1", "k2", "k3",
                                                      "k4", "k5", "k6", "k7","k8","k9",
                                                      "k10"))
  return(new_dataset)
}
# ----- [1.3] dataset importation ----------------------------------------------

# Import the data from a tab delimited ascii file 
data<-read.table(file= 'gamma_selection.txt',
                 header=T,
                 dec=".", fill = T) 
names(data)
data$rate=ordered(data$rate, levels=c("k1", "k2", "k3",
                                    "k4", "k5", "k6", "k7","k8","k9",
                                    "k10"))
colnames(data) = c("rate","likelihood", "lambda", "alpha")
summary(data)
data = subset(data,data$lambda!="NA")

for (u in c(1:10)){ 
  ku = paste0("k",u)
  data$AIC[data$rate==ku]= 2*u - 2*(-1*data$likelihood[data$rate==ku]) 
}




# ----- [1.3] lambda plot - general --------------------------------------------

ggplot(data = data, aes(rate,lambda))+
  geom_boxplot(aes(fill=rate),  outlier.shape = NA, show.legend = F )+
  geom_jitter(alpha=0.3, height = 0)

ggsave(filename = "gamma_selection_eusoc_lambda.tiff",
       width = 7,height = 8,units = "in", 
       dpi = 320, limitsize = F)

ggplot(data = data, aes(rate,AIC))+
  stat_boxplot(aes(fill=rate) ,outlier.shape = NA,show.legend = F )+
  geom_jitter(alpha=0.3)+
  ylab("AIC")

ggsave(filename = "gamma_selection_plot_AIC.tiff",
       width = 7,height = 8,units = "in", 
       dpi = 320, limitsize = F)


# ----- [1.3] lambda plot - only converging ------------------------------------


new_data = keep_median_l_values(data)

ggplot(data = new_data, aes(rate,lambda))+
  geom_boxplot(aes(fill=rate),  outlier.shape = NA, show.legend = F )+
  geom_jitter(alpha=0.3)


ggsave(filename = "gamma_selection_plot_med_lambda.tiff",
       width = 7,height = 8,units = "in", 
       dpi = 320, limitsize = F)

ggplot(data = new_data, aes(rate,-likelihood))+
  stat_boxplot(aes(fill=rate) ,outlier.shape = NA,show.legend = F )+
  geom_jitter(alpha=0.3)+
  ylab("ln(likelihood)")

ggsave(filename = "gamma_selection_plot_med_lnL.tiff",
       width = 7,height = 8,units = "in", 
       dpi = 320, limitsize = F)


# Looking at the graph, we will keep k=1,2,3. This is Because:
# - lambda values are converging for k = 1 and k = 2
# - the lowest AIC is for k=3 so we can test it but it may not converge

mean(new_data$rate[new_data$rate=="k2"])

# ln(likelihood)[k=2] = -325,959
mean(new_data$likelihood[new_data$rate=="k1"])
# ln(likelihood)[k=1] = -380,918