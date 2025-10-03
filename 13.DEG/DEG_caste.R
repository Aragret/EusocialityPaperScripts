################################################################################
#
#             Title: Differential gene expression 
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
The aim of theanalysis is to compare the gene expresion of each caste in a 
pairwise manner to identify genes that are caste specific. The output should
be an excel/csv table containing the gene names, the bias for each pairwise 
comparison, and the overall specificty.
e.g.  gene_name Wo_So Wo_Al So_Al Specificity
      OG.000001 Wo    Wo    So    none
      OG.000002 Wo    Wo    none  Wo
      OG.000003 Wo    none  Al    none
'


# ----- [0.2] Dependancies: ----------------------------------------------------

'The following libraries are required: DESeq2,  dplyr, ggplot2'

# ------------------------------------------------------------------------------
# ------ [1] Introduction ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- [1.1] Load packages ----------------------------------------------------

library(DESeq2)
library(dplyr)
library(ggplot2)

# ----- [1.2] Set directory ----------------------------------------------------

setwd("C:/Users/caumont/Documents/WORK/Dino_McMahon/Work/Lab work and projects/6. Gene network/6. DEG/1.DESEQ2/script/")

# ----- [1.3] prepare output ----------------------------------------------------

# chose the type of dataset to process
input_reads = "beforeBC" # or "beforeBC" afterBC
type_output = "LFC" # or "LFCval"

df_out = data.frame(Species = character(),
  gene_name = character(),
	Specificity = character(),
	Al_Ne = character(),
	Al_Ny = character(),
	Al_So = character(),
	Al_Wo = character(),
	Ne_Ny = character(),
	Ne_So = character(),
	Ne_Wo = character(),
	Ny_So = character(),
	Ny_Wo = character(),
	So_Wo = character(),
	Er_So = character(),
	Er_Wo = character(),
	Nm_Ny = character(),
	Nm_So = character(),
	Nm_Wo = character(),
	Pr_So = character(),
	Pr_Wo = character(),
	Ad_Ju = character()
)



# ------------------------------------------------------------------------------
# ------ [2] Main --------------------------------------------------------------
# ------------------------------------------------------------------------------


# ------ Preparation of lists --------------------------------------------------
# Species list
species_list = c(               # Attention: the order matters: need to be the 
                                #same as spe_list
  "Anoplotermes_pacificus",
  "Coptotermes_gestroi",
  "Hodotermopsis_sjostedti",
  "Kalotermes_flavicollis",
  "Mastotermes_darwiniensis",
  "Macrotermes_natalensis",
  "Neotermes_castaneus",
  "Prorhinotermes_simplex",
  "Reticulitermes_flavipes",
  "Zootermopsis_nevadensis",
  "Cryptocercus_meridianus",
  "Cryptocercus_punctulatus",
  "Blattella_germanica",
  "Blatta_orientalis")

# Abbreviation list
spe_list = c("Apac","Cges","Hsjo","Kfla", "Mdar","Mnat","Ncas","PRsim","Rfla",
             "Znev","Cmer","Cpun","Bger","Bori")



for (SPECIES in species_list){
  
  spe = spe_list[match(SPECIES,species_list)]




# ----- [1.] Upload read count ------------------------------------------------
  readcounts_file = paste0("../data/normalised_count_",
                           input_reads,
                           "/",
                           SPECIES,
                           ".nreadcount.txt")
  readcounts=read.table(readcounts_file,
                        header=TRUE,
                        row.names=1,
                        sep="\t",
                        stringsAsFactors=FALSE)
  dim(readcounts)

# ----- [1.] Upload conditions -------------------------------------------------
  conditions_file = paste0("../data/conditions/",
                           SPECIES,
                           ".conditions.txt")
  coldata=read.table(conditions_file,
                     header=TRUE,
                     row.names=1, 
                     sep="\t", 
                     stringsAsFactors=FALSE)
  dim(coldata)

  readcounts = readcounts %>% relocate(rownames(coldata))

# -----  factors ---------------------------------------------------------------

  coldata$caste = factor(coldata$caste)
  coldata$sex = factor(coldata$sex)
  
# ------------------------------------------------------------------------------
# ------ [2] DEG ---------------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- [2.2] Model ------------------------------------------------------------

  dds=DESeqDataSetFromMatrix(countData = readcounts,
                           colData = coldata,
                           design = ~ caste)

# ----- [2.3] pre-filtering ----------------------------------------------------

# Filtering the very low expressed genes
  keep = rowSums(counts(dds)) >= 0
  dds = dds[keep,]
  dds

# ----- [2.4] DESeq object -----------------------------------------------------

  dds = DESeq(dds)
  
# ------- choice of the caste to compare ---------------------------------------

  Fin=data.frame(row.names(readcounts))
  colnames(Fin) = "gene_name" 
  headers=c("gene_name")
  
#loop here, for casteA, casteB in level coldata$caste...
  for (i in c(1:(length(levels(coldata$caste))-1))) {
    casteA = levels(coldata$caste)[i]
    for (j in c((i+1):length(levels(coldata$caste)))){
      casteB = levels(coldata$caste)[j]
      print(casteA)
      print(casteB)
  
      # dds result computation
      result = results(dds,contrast = c("caste", casteA, casteB))
      df = data.frame(result)
      df[is.na(df)] = 1
      setpval = 0.05
      df$DEG = "none"
      
      # changing the Er into Ne for 3 species to save the correct caste in Fin df
      if (spe %in% c("Kfla", "Ncas", "Hsjo")){
        if (casteA == "Er"){casteA = "Ne"}
        if (casteB == "Er"){casteB = "Ne"}
      }
      
      for (line in c(1:dim(df)[1])) {
        if (df[line,6] < setpval & df[line,2] > 0) {
          if (type_output == "LFC"){df$DEG[line]=casteA}
          if (type_output == "LFCval"){df$DEG[line]=round(df[line,2],2)}
        }
        
        if (df[line,6] < setpval & df[line,2] < 0) {
          if (type_output == "LFC"){df$DEG[line]=casteB}
          if (type_output == "LFCval"){df$DEG[line]=round(df[line,2],2)}
        }
        
      }
      colnames(df)[dim(df)[2]]=paste0(casteA, "_",casteB)
      headers = append(headers,paste0(casteA, "_",casteB))
      # Saving in Fin df
      Fin = cbind.data.frame(Fin, df[,dim(df)[2]])

      # replacing back the caste for the dds result computation
      if (spe %in% c("Kfla", "Ncas", "Hsjo")){
        if (casteA == "Ne"){casteA = "Er"}
        if (casteB == "Ne"){casteB = "Er"}
      }
      
    }
  }

  head(Fin)
  colnames(Fin)=headers

  if (type_output == "LFC"){
    Fin$Specificity="NA" 
    goal = (length(levels(coldata$caste))-1)
    all = dim(Fin)[2]-2
    for (i in c(1:dim(Fin)[1])){
      Al = 0
      Ny = 0
      Er = 0
      Pr = 0
      Wo = 0
      So = 0
      Nm = 0
      Ad = 0
      Ju = 0
      Ne = 0
      none = 0
      for (j in Fin[i,]) {
        if (j=="Al") { Al = Al+1}
        if (j=="Ny") { Ny = Ny+1}
        if (j=="Er") { Er = Er+1}
        if (j=="Pr") { Pr = Pr+1}
        if (j=="Wo") { Wo = Wo+1}
        if (j=="So") { So = So+1}
        if (j=="Ju") { Ju = Ju+1}
        if (j=="Ad") { Ad = Ad+1}
        if (j=="Nm") { Nm = Nm+1}
        if (j=="Ne") { Ne = Ne+1}
        if (j=="none") { none = none+1}
      }
      if (none == (all - goal)){
        if (Al == goal){Fin$Specificity[i]="Al"}
        if (Ny == goal){Fin$Specificity[i]="Ny"}
        if (Er == goal){Fin$Specificity[i]="Er"}
        if (Pr == goal){Fin$Specificity[i]="Pr"}
        if (Wo == goal){Fin$Specificity[i]="Wo"}
        if (So == goal){Fin$Specificity[i]="So"}
        if (Ju == goal){Fin$Specificity[i]="Ju"}
        if (Ad == goal){Fin$Specificity[i]="Ad"}
        if (Nm == goal){Fin$Specificity[i]="Nm"}
        if (Ne == goal){Fin$Specificity[i]="Ne"}
      }
    }
  }
  Fin$Species = spe
  df_out = bind_rows(df_out, Fin)
}

# ------------------------------------------------------------------------------
# ------ Write results ---------------------------------------------------------
# ------------------------------------------------------------------------------


write.csv2(df_out,
           paste0("../result/DEG_",
                  input_reads,
                  "_per_caste_only_",
                  type_output,
                  "_n005_FINAL.csv"))

