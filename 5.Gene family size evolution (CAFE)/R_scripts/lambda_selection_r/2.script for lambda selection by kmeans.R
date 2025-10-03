################################################################################
#
#             Title: Lambda k-means clustering 
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
This script is used to cluster the evolutionnary rates "lambda" calculated by 
CAFE5 each branch of the termite phylogenetic tree.

In this script, we use tree tests to identify the optimal number of cluster 
for lambda "km". After running the tree tests (elbow, silhouette, gap stat),
we display the PCA including the background and the foreground lambda.
'


# ----- [0.2] Dependancies: ----------------------------------------------------

'The following libraries are required: factoextra, FactoMineR, ggplot2, cluster'
'The following input files are required: Mlambda_selection.txt
The output files are:
clusters2.csv
clusters3.csv
clusters4.csv
clusters5.csv
clusters6.csv
'
# ------------------------------------------------------------------------------
# ------ [1] Introduction ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- [1.1] Load packages ----------------------------------------------------

library(factoextra)
library(FactoMineR)
library(ggplot2)
library(cluster)

# ----- [1.2] Set directory --------------------------------------------------------------

# Please set to data path on your system
pa = "C:/Users/caumont/Documents/WORK/Dino_McMahon/Work/Lab work and projects/"
setwd(
  paste0(
    pa, "10. Eusociality/CAFE/4. Selection"
  )
)

# ----- [1.3] Load data ------------------------------------------------------------------


data <- read.table('Mlambda_selection.txt', #Mk2_lambda_selection.txt
                   dec = ".",
                   header = T,
                   stringsAsFactors = T,
                   fill = T)
str(data)#row.names(data) = data$Branch


lambda2 = data[,c(4)]
boxplot(lb2~mod, data = data)
colnames(data) = c("Branch", "lnL", "lambda1", "lambda2")
data = subset(data,data$lambda2!="NA")
# ------------------------------------------------------------------------------
# ------ [2] lambda convergence ------------------------------------------------
# ------------------------------------------------------------------------------


convergence_test = function(data, test){ 
  lambda0_list = c()
  lambda_list = c()
  M_list = c()
  not_conv = c()
  
  for (M in levels(data$Branch)){
    B_M1_values = data[data$Branch %in% c(M),4]
    B_M1_rounds = round(B_M1_values,3)
    B_M1_rndmed = median(B_M1_rounds)
    B_M1_median = median(B_M1_values)
    B_M1_select = B_M1_values[B_M1_rounds==B_M1_rndmed]
    half_len = length(B_M1_values)*0.5
    len_select = length(B_M1_select)
    if (test){
      if (len_select > half_len){
        print(paste0("Branch ", M," converged"))
        lambda0_list = append(lambda0_list, median(data[data$Branch %in% c(M),3]))
        lambda_list = append(lambda_list, B_M1_median)
        M_list = append(M_list, M)
      }
      else{
        print(paste0("Branch ", M," did not converged"))
        not_conv = append(not_conv, M)
      }
    }
    else{
      info = "not a test"
      lambda0_list = append(lambda0_list, median(data[data$Branch %in% c(M),3]))
      lambda_list = append(lambda_list, B_M1_median)
      M_list = append(M_list, M)
    }
  }
  new_dataf = data.frame(cbind(M_list,lambda0_list,lambda_list)) 
  colnames(new_dataf)= c("Branch", "lambda0","lambda")
  if (test){return(not_conv)}
  if (!test){print(info)
    return(new_dataf)}
}

# ------------------------------------------------------------------------------

M_not_conv =convergence_test(data, T)

# Check the branches that did not converge and choose what to do
ggplot(data = data[data$Branch %in% c("M50"),],
       aes(Branch, lambda2)) +
  geom_boxplot(aes( fill=Branch), 
               outlier.shape=NA) +
  geom_jitter() +
  ggtitle(label = "Branches that did not converge",
          subtitle = "What values should we pick?")

# Here we decided to take the median of the values

med_df = convergence_test(data, F)

row.names(med_df) = med_df$Branch
med_df$lambda0=as.numeric(med_df$lambda0)
med_df$lambda=as.numeric(med_df$lambda)

conv_lambda = as.numeric(med_df[,c(3)])

# ------------------------------------------------------------------------------
# ------ [2] Finding km --------------------------------------------------------
# ------------------------------------------------------------------------------


# ----- [2.1] Elbow test function ----------------------------------------------


elbow = function(data){

  set.seed(123)
  wss_values=c()
  k.values <- 1:15        # Compute and plot wss for k = 1 to k = 15
  
  for (i in k.values){ # extract wss for 2-15 clusters
    wss = kmeans(data, i, nstart = 1000)$tot.withinss # compute total within-cluster sum of square
    wss_values <- append(wss_values, wss)
  }
  
  plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

}


# ----- [2.2] Silhouette test --------------------------------------------------

Silhouette = function(data){
  avg_sil_values <- c()
  k.values <- 2:15
  for (i in k.values){
    res <- kmeans(data, centers = i, nstart = 1000)
    ss <- silhouette(res$cluster, dist(data))
    avg_sil <- mean(ss[, 3])
    avg_sil_values <- append(avg_sil_values, avg_sil)
  }
  plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")
}



# ----- [2.3] gap stat ---------------------------------------------------------

gap_stat = function(data){
  g_stat <- clusGap(data[,c(3,3)], FUN = kmeans, nstart = 1000,
                    K.max = 15, B = 50)
  # Print the result
  fviz_gap_stat(g_stat)
}


# ---- [2.4] Tests -------------------------------------------------------------

elbow(conv_lambda)
Silhouette(conv_lambda)
plot(sort(conv_lambda))

# ------------------------------------------------------------------------------
# ------ [3] Clustering --------------------------------------------------------
# ------------------------------------------------------------------------------


clustering = function(data, km){
res = kmeans(data,
             centers = km,
             iter.max = 100,
             nstart =  1000,
             algorithm =  "Hartigan-Wong",
             trace =  FALSE)

fviz_cluster(res, data = med_df[,c(3,3)])
}
clustering(med_df[,c(3,3)], 2)
clustering(med_df[,c(3,3)], 3)
clustering(med_df[,c(3,3)], 4)
clustering(med_df[,c(3,3)], 5)
clustering(med_df[,c(3,3)], 6)


# if lambda 1 and 2
res.pca = PCA(res$cluster)
plot.PCA(res.pca)

# write cluster for optimal tree
res = kmeans(med_df[,c(3,3)],
             centers = 6,
             iter.max = 100,
             nstart =  1000,
             algorithm =  "Hartigan-Wong",
             trace =  FALSE)
fviz_cluster(res, data = med_df[,c(3,3)])
write.csv2(res$cluster, "clusters.csv")
df=read.csv2("clusters.csv")
colnames(df)= c("Branch","Cluster")
for (i in c(1:9)){
  M_ = paste0("M",i)
  M0 = paste0("M0",i)
  df$Branch[df$Branch==M_] = M0
}
write.csv2(df, "clusters5.csv")




