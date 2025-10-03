#
#   Title: Glmm on CAFE count
#
#
#############################################################

#-------------------------------------------------------------------------------
#---------------------------- Libraries -------------------------------
#-------------------------------------------------------------------------------



library(tidyverse)
library(ape)
library(geiger)
library(MCMCglmm)
library(dplyr)
# library(coda) 

setwd("C:/Users/caumont/Documents/WORK/Dino_McMahon/Work/Lab work and projects/10. Eusociality/MCMCglmm_CAFE")


#-------------------------------------------------------------------------------
#---------------------------- Preparing datasets -------------------------------
#-------------------------------------------------------------------------------

data = read.csv("data/dataset_mcmcglmm_cafe.csv",
                header = T,
                sep = ";",
                dec = ",", 
                stringsAsFactors = T)

str(data)




# ------------ checking the HOGs ---------------------------------------------


fin_df = data.frame(hog = character(),
                    p_val = numeric(),
                    p_adj = numeric(),
                    covariate = character(),
                    prior = character(),
                    randon = character())
for (var1 in c("Ontogeny", "Social_category")){
  for (prior in c("iw", "expanded")){
    for (rd in c("with", "without")){
      depvar=paste0(var1," with ",prior," priors and ",rd," random")
      filenames = list.files("result_folder", pattern = paste0("*",depvar,"*"), full.names = T) 
      N_row = length(filenames)
      j=0
      p_list=c()
      hog_list=c()
      for (modfile in filenames){
        print(modfile)
        j = j + 1
        mod = load(modfile)
        model = eval(as.name(mod))
        hog_name = strsplit(modfile, " ")[[1]][5]
        hog_list=append(hog_list, hog_name)
        p_temp = summary(model)$solutions[2,5]
        p_list=append(p_list, p_temp)
      }
      p_adj = p.adjust(p_list,method = "fdr")
      df=data.frame(list(hog = hog_list,
                         p_val = p_list,
                         p_adj = p_adj,
                         covariate = var1,
                         prior = prior,
                         randon = rd))
      df_psig = df[df$p_adj<0.05,]
      dim(df_psig)
      fin_df = bind_rows(fin_df,df_psig)
    }}}

write.csv2(x = fin_df, file = paste0("Result_Disc/Significance sigHOG.csv"),
           quote = F, row.names = F)


# ----------- Quick figures ----------------------------------------------------


# Ontogeny

depvar="Ontogeny with expanded priors and with random"
filenames = list.files("result", pattern = paste0("*",depvar,"*"), full.names = T) 
N_row = length(filenames)
j=0
p_list=c()
hog_list=c()
for (modfile in filenames){
  print(modfile)
  j = j + 1
  mod = load(modfile)
  model = eval(as.name(mod))
  hog_name = strsplit(modfile, " ")[[1]][5]
  hog_list=append(hog_list, hog_name)
  p_temp = summary(model)$solutions[2,5]
  p_list=append(p_list, p_temp)
}
p_adj = p.adjust(p_list,method = "fdr")
df=data.frame(list(hog = hog_list,
                   p_val = p_list,
                   p_adj = p_adj))
df_psig = df[df$p_adj<0.05,]
dim(df_psig)



pdf("Result_Disc/sigHOG ontogeny  - quick figure.pdf",width = 7,height = 9)
par(mfrow=c(3,3))
for (no_hog in unique(fin_df$hog[fin_df$covariate=="Ontogeny"])){
  #no_hog = "N0.HOG0000380"
  boxplot(data[data$Family_ID==no_hog,]$gene_count~data[data$Family_ID==no_hog,]$Ontogeny, xlab = no_hog)
  points(data[data$Family_ID==no_hog,]$gene_count~data[data$Family_ID==no_hog,]$Ontogeny, col="blue",pch=19)
  points(data[data$Family_ID==no_hog & data$Species=="Mdar",]$gene_count~data[data$Family_ID==no_hog & data$Species=="Mdar",]$Ontogeny,
         pch=19, col = "red")
  points(data[data$Family_ID==no_hog & data$Species=="Dlon",]$gene_count~data[data$Family_ID==no_hog & data$Species=="Dlon",]$Ontogeny,
         pch=19, col = "green")
}
dev.off()


# Social_category

depvar="Social_category with iw priors and with random"
filenames = list.files("result", pattern = paste0("*",depvar,"*"), full.names = T) 
N_row = length(filenames)
j=0
p_list=c()
hog_list=c()
for (modfile in filenames){
  print(modfile)
  j = j + 1
  mod = load(modfile)
  model = eval(as.name(mod))
  hog_name = strsplit(modfile, " ")[[1]][5]
  hog_list=append(hog_list, hog_name)
  p_temp = summary(model)$solutions[2,5]
  p_list=append(p_list, p_temp)
}
p_adj = p.adjust(p_list,method = "fdr")
df=data.frame(list(hog = hog_list,
                   p_val = p_list,
                   p_adj = p_adj))
df_psig = df[df$p_adj<0.05,]
dim(df_psig)


pdf("Result_Disc/sigHOG Social_category  - quick figure.pdf",
    width = 7,height = 9)
par(mfrow=c(3,3))
for (no_hog in unique(fin_df$hog[fin_df$covariate=="Social_category"])){
  #no_hog = "HOG0000000"
  boxplot(data[data$Family_ID==no_hog,]$gene_count~data[data$Family_ID==no_hog,]$Social_category, xlab = no_hog)
  points(data[data$Family_ID==no_hog,]$gene_count~data[data$Family_ID==no_hog,]$Social_category, col="blue",pch=19)
  points(data[data$Family_ID==no_hog & data$Species=="Mdar",]$gene_count~data[data$Family_ID==no_hog & data$Species=="Mdar",]$Social_category,
         pch=19, col = "red")
  points(data[data$Family_ID==no_hog & data$Species=="Dlon",]$gene_count~data[data$Family_ID==no_hog & data$Species=="Dlon",]$Social_category,
         pch=19, col = "green")
}
dev.off()



for (no_hog in unique(Hog)){
#no_hog = "N0.HOG0000380"
boxplot(data[data$Family_ID==no_hog,]$gene_count~data[data$Family_ID==no_hog,]$Ontogeny, xlab = no_hog)
points(data[data$Family_ID==no_hog,]$gene_count~data[data$Family_ID==no_hog,]$Ontogeny, col="blue",pch=19)
points(data[data$Family_ID==no_hog & data$Species=="Mdar",]$gene_count~data[data$Family_ID==no_hog & data$Species=="Mdar",]$Ontogeny,
       pch=19, col = "red")
points(data[data$Family_ID==no_hog & data$Species=="Dlon",]$gene_count~data[data$Family_ID==no_hog & data$Species=="Dlon",]$Ontogeny,
       pch=19, col = "green")
}
write.csv2(x = df_psig, file = "hogsize_following_ontogeny_pattern.csv",
          quote = F, row.names = F)
Hog = c("N0.HOG0000380","N0.HOG0000481","N0.HOG0000502","N0.HOG0000906",
        "N0.HOG0000937","N0.HOG0000964","N0.HOG0000974","N0.HOG0000976",
        "N0.HOG0001119","N0.HOG0001245","N0.HOG0001248","N0.HOG0001292",
        "N0.HOG0001513","N0.HOG0001526","N0.HOG0001575","N0.HOG0001671",
        "N0.HOG0001738","N0.HOG0001798","N0.HOG0001834","N0.HOG0001851",
        "N0.HOG0001878","N0.HOG0002107","N0.HOG0002353","N0.HOG0002845",
        "N0.HOG0003077","N0.HOG0003410","N0.HOG0000380","N0.HOG0000481",
        "N0.HOG0000502","N0.HOG0000659","N0.HOG0000781","N0.HOG0000906",
        "N0.HOG0000937","N0.HOG0000964","N0.HOG0000976","N0.HOG0001119",
        "N0.HOG0001245","N0.HOG0001248","N0.HOG0001292","N0.HOG0001513",
        "N0.HOG0001526","N0.HOG0001537","N0.HOG0001575","N0.HOG0001671",
        "N0.HOG0001738","N0.HOG0001798","N0.HOG0001834","N0.HOG0001851",
        "N0.HOG0001878","N0.HOG0001979","N0.HOG0002048","N0.HOG0002107",
        "N0.HOG0002323","N0.HOG0002353","N0.HOG0002845","N0.HOG0002865",
        "N0.HOG0002943","N0.HOG0003014","N0.HOG0003077","N0.HOG0003410",
        "N0.HOG0004986","N0.HOG0006820","N0.HOG0000380","N0.HOG0000481",
        "N0.HOG0000502","N0.HOG0000659","N0.HOG0000906","N0.HOG0000937",
        "N0.HOG0000964","N0.HOG0000974","N0.HOG0000976","N0.HOG0001245",
        "N0.HOG0001248","N0.HOG0001292","N0.HOG0001513","N0.HOG0001526",
        "N0.HOG0001671","N0.HOG0001738","N0.HOG0001798","N0.HOG0001834",
        "N0.HOG0001851","N0.HOG0001878","N0.HOG0001979","N0.HOG0002048",
        "N0.HOG0002107","N0.HOG0002353","N0.HOG0002845","N0.HOG0003014",
        "N0.HOG0003077","N0.HOG0003410","N0.HOG0000380","N0.HOG0000481",
        "N0.HOG0000502","N0.HOG0000659","N0.HOG0000906","N0.HOG0000937",
        "N0.HOG0000964","N0.HOG0000974","N0.HOG0000976","N0.HOG0001119",
        "N0.HOG0001245","N0.HOG0001248","N0.HOG0001292","N0.HOG0001513",
        "N0.HOG0001526","N0.HOG0001575","N0.HOG0001671","N0.HOG0001738",
        "N0.HOG0001798","N0.HOG0001834","N0.HOG0001851","N0.HOG0001878",
        "N0.HOG0001979","N0.HOG0002048","N0.HOG0002107","N0.HOG0002353",
        "N0.HOG0002845","N0.HOG0003014","N0.HOG0003077","N0.HOG0003410",
        "N0.HOG0000380","N0.HOG0000481","N0.HOG0000502","N0.HOG0000659",
        "N0.HOG0000906","N0.HOG0000937","N0.HOG0000964","N0.HOG0000974",
        "N0.HOG0000976","N0.HOG0001119","N0.HOG0001245","N0.HOG0001248",
        "N0.HOG0001292","N0.HOG0001476","N0.HOG0001513","N0.HOG0001526",
        "N0.HOG0001575","N0.HOG0001671","N0.HOG0001738","N0.HOG0001798",
        "N0.HOG0001834","N0.HOG0001851","N0.HOG0001878","N0.HOG0001979",
        "N0.HOG0002048","N0.HOG0002107","N0.HOG0002350","N0.HOG0002353",
        "N0.HOG0002845","N0.HOG0002943","N0.HOG0003014","N0.HOG0003077",
        "N0.HOG0003085","N0.HOG0003410","N0.HOG0004986")
unique(Hog)
length(unique(Hog))


interesting_hogs = c("N0.HOG0000380", "N0.HOG0000937", "N0.HOG0001513",
                     "N0.HOG0001851", "N0.HOG0003410", "N0.HOG0002865",
                     "N0.HOG0002943", "N0.HOG0004986", "N0.HOG0006820")

data$Ontogeny= factor(data$Ontogeny,
                      levels = c("linear", "bifurcated"),
                      labels = c("Species with\nlinear ontogeny",
                                 "Species with\nbifurcated ontogeny"))
data$specol = NA
data$specol[data$Species=="Mdar"] = "M. darwiniensis"
data$specol[data$Species=="Dlon"] = "D. longilabius"
data$specol[data$Species=="Cmer"] = "C. meridianus"
data$specol[data$Species=="Cpun"] = "C. punctulatus"

p <- ggplot(data = data[data$Family_ID %in% interesting_hogs,],
            aes(x=Ontogeny, y=gene_count)) 
p <- p + geom_boxplot(aes(fill=Ontogeny), outlier.colour = NA)
p <- p + geom_jitter(aes(shape = specol), width = 0.1, height = 0, size = 1.5)
p <- p + facet_wrap( ~ Family_ID, scale="free")
p <- p + ylab("Gene count") #+ ggtitle("Title")
p <- p + xlab("")
p <- p + guides(fill=guide_legend(title="Ontogeny",keyheight = 2,title.hjust = .3))
# p <- p + scale_shape_discrete(solid = T, na.translate = F)
p <- p + scale_shape_manual(values=c(15,16,17,18), na.translate = F,
                            guide= guide_legend(label.theme = element_text(face = "italic",size = 10),
                                                title = "Species", title.hjust = .3))
p <- p + scale_x_discrete(labels = c('',''))
p <- p + scale_fill_manual( values = c("mediumpurple3", "steelblue1"))
#p <- p + scale_fill_brewer(palette="OrRd")
p <- p + theme(axis.ticks = element_blank(), 
               panel.background = element_rect(fill = 'white', colour = 'white'),
               panel.grid.minor.y = element_line(colour = "grey80",linetype = "dotted"),
               panel.grid.major.y = element_line(colour = "grey",linetype = "dotted"),
               panel.grid.major.x =  element_line(colour = "grey",linetype = "dotted"),
               plot.margin = unit(c(.5, 0, .5, 0), "cm"),
               axis.line.x  = element_line(color = "black"),
               axis.line.y  = element_line(color = "black"),
               axis.title.x = element_text(vjust = -1)) 
p
ggsave("Figure_Selected_hogs_pattern_ontogeny.png",width = 7.5,height = 7, dpi=660)


# ------------ getteing HOG accepted function ----------------------------------

HOG_accepted = c("HOG0000380",
                 "HOG0001851",
                 "HOG0002046",
                 "HOG0002649",
                 "HOG0004986",
                 "HOG0000937",
                 "HOG0001513",
                 "HOG0002943",
                 "HOG0006820",
                 "HOG0000758")

library(readxl)
master = data.frame(read_excel("../../6. Gene network/TABLES/Master_sheet_eusoc.xlsx",
                               sheet = "S0_Master_sheet",
                               na = c("","NA")))

info_HOGacc = subset(master, subset = hog %in% HOG_accepted)

write.csv2(info_HOGacc, "Result_Disc/information_sigHOG_following_sociality.csv",
           row.names = F,quote = F)

