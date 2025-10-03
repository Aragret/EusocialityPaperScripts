################################################################################

# Author: Cedric Aumont
# Year: 2025

################################################################################


#-------------------------------------------------------------------------------
#-------Introduction -----------------------------------------------------------
#-------------------------------------------------------------------------------

"%ni%" <- Negate("%in%")
len = length
library(wTO)
library(readxl)
library(dplyr)


#set wd to source file location
#setwd("~/WORK/Dino_McMahon/Work/Lab work and projects/6. Gene network/OR_to_DEG")


# upload data



spe_list = c("Apac","Cges","Hsjo","Kfla", "Mdar","Mnat","Ncas","PRsim","Rfla")

species_list = c(               # Attention: the order matters: need to be the same as spe_list
  "Anoplotermes_pacificus",
  "Coptotermes_gestroi",
  "Hodotermopsis_sjostedti",
  "Kalotermes_flavicollis",
  "Mastotermes_darwiniensis",
  "Macrotermes_natalensis",
  "Neotermes_castaneus",
  "Prorhinotermes_simplex",
  "Reticulitermes_flavipes")


ORs112_g = read.table("../Data/ORs112g_spe_gene_caste_correct.txt", sep="\t",header = T)
allOG_g = read.table("../Data/allOGg_spe_gene_caste_correct.txt",header = T, sep="\t")
clus2_g = read.table("../Data/clus2g_spe_gene_caste_correct.txt",header = T, sep="\t")

alldeg_wosore = read.table("../Data/DEGwosore_genelist_termites.txt", 
                           col.names = c("spe", "gene"))

sHOG_DEGpWORE = read_xlsx(path = "../Data/Results_pDEG_WORE.xlsx",
                          sheet = "sHOG_DEGpWORE-nostat")

TW = read.table("../Data/DEG_TW.txt", sep="\t",header = T)$x
FW = read.table("../Data/DEG_FW.txt", sep="\t",header = T)$x
TR = read.table("../Data/DEG_TR.txt", sep="\t",header = T)$x
FR = read.table("../Data/DEG_FR.txt", sep="\t",header = T)$x


#-------------------------------------------------------------------------------
#------- Figure -----------------------------------------------------
#-------------------------------------------------------------------------------




#Figure

OR112orall = "all"

for (spe in spe_list){  
  for (casteNtw in c("Wo", "Re")){
    for (casteOr in c("Wo", "Re")){
      #for (casteDEG in c("Wo", "Re")){
      # choose title network
      title = paste0("../Figures_output/output/Really_Corrected_Figure_5_",casteNtw,"NtW_",casteOr,"OR_","wore_DEG_",spe, "_",OR112orall)
      titletxt = paste0("Really_Corrected_Figure_5_",casteNtw,"NtW_",casteOr,"OR_","wore_DEG_",spe, "_",OR112orall,".txt")
      species = species_list[match(spe,spe_list)]
      if(OR112orall == "all"){
        network = paste0("../Networks/oralldeg_", species, "_",casteNtw,"_nn1.RData")
        load(network)
        if (casteOr == "Wo"){
          ORs = allOG_g$gene[allOG_g$spe == spe & allOG_g$caste == "Wo"]
        } else {
          ORs = allOG_g$gene[allOG_g$spe == spe & allOG_g$caste == "Re"]}
      }
      # choose the DEG type
      DEG_spe = sHOG_DEGpWORE[sHOG_DEGpWORE$spe == spe,]
      onto_deg =c(TW, TR,FW, FR)

      DEG_bia = intersect(subset(DEG_spe, gene %in% onto_deg)$gene, alldeg_wosore$gene)
      
      
      Ntw_OD=
        Network[Network$Node.1 %in% c(ORs, DEG_bia) & Network$Node.2 %in% c(ORs, DEG_bia),]
      
      checkntw = Ntw_OD[abs(Ntw_OD$wTO)>=.4,] 
      in_ntw = unique(c(checkntw$Node.1,checkntw$Node.2))
    #  non_OG112 = SuperExactTest::intersect(
    #    allOG$gene[allOG$spe == spe & allOG$gene %ni% OG112$gene],
    #    in_ntw)
    #  my_OG112 = SuperExactTest::intersect(
    #    OG112$gene[OG112$spe == spe],
    #    in_ntw)                                   
      if (dim(checkntw)[1] >2){
        #Vis = NetVis2(Node.1 = Ntw_OD$Node.1,
         #             Node.2 = Ntw_OD$Node.2, pval = NULL,
          #            wTO = Ntw_OD$wTO,MakeGroups ="manual",#walktrap
           #           Cluster =F,cutoff = list(kind="Threshold", value = 0.4),
            #          smooth.edges = F,
     #                 shape = list(shape=c(rep("triangle",len(my_OG112)),rep("square",len(non_OG112))),
    #                             names = c(my_OG112, non_OG112)),
                   #   layout = "layout_with_kk")#,
        #path = paste0("Network_ORtoDEG2/",title ,".html"))
        write.table(Ntw_OD, file = paste0("../Figures_output/output/",titletxt),quote = F,row.names = F,dec = ".", sep="\t")
        
      }
    }
  }
}




#nodes$shape =  "dot"



# ------- Netvis function modified ---------------------------------------------

NetVis2=function (Node.1, Node.2, wTO, pval = NULL, MakeGroups = FALSE, 
                  padj = NULL, cutoff = list(kind = "Threshold", value = 0.5), 
                  layout = NULL, smooth.edges = T, path = NULL, Cluster = F, 
                  legend = T, shape = list(shape = "triangle", names = NULL), 
                  manipulation = F) 
{
  input_vis = data.frame(Node.1 = Node.1, Node.2 = Node.2, 
                         wTO = as.numeric(wTO))
  if (!is.null(pval)) {
    input_vis$pval = pval
  }
  if (!is.null(padj)) {
    input_vis$padj = padj
  }
  `%ni%` <- Negate(`%in%`)
  `%>%` <- magrittr::`%>%`
  if (cutoff$kind %ni% c("Threshold", "pval", "pval.adj")) {
    stop("cutoff kind must be \"Threshold\", \"pval\" or \"pval.adj\".")
  }
  if (is.numeric(cutoff$value) == F) {
    stop("cutoff value must be numeric.")
  }
  MakeGroups_pos = c("walktrap", "optimal", "spinglass", "edge.betweenness", 
                     "fast_greedy", "infomap", "louvain", "label_prop", "leading_eigen", 
                     FALSE, "manual")
  if (MakeGroups %ni% MakeGroups_pos) {
    stop("MakeGroups should be FALSE or one of the following options: 'walktrap',\n           'optimal', 'spinglass', 'edge.betweenness',\n           'fast_greedy', 'infomap', 'louvain', 'label_prop',\n           'leading_eigen'.")
  }
  if (Cluster %ni% c(T, F)) {
    stop("Cluster must be T / F.")
  }
  if (smooth.edges %ni% c(T, F)) {
    stop("smooth.edges must be T / F.")
  }
  input_vis = subset(input_vis, abs(input_vis$wTO) > 0.01)
  if (cutoff$kind == "Threshold") {
    input_vis = subset(input_vis, abs(input_vis$wTO) >= 
                         cutoff$value)
  }
  else if (cutoff$kind == "pval") {
    input_vis = subset(input_vis, input_vis$pval <= cutoff$value)
  }
  else if (cutoff$kind == "pval.adj") {
    input_vis = subset(input_vis, input_vis$pval.adj <= 
                         cutoff$value)
  }
  if (nrow(input_vis) <= 2) {
    stop("There is less than 2 nodes on your network. Choose a lower cutoff.")
  }
  if (smooth.edges == T) {
    smooth.edges = "enabled"
  }
  input_vis = input_vis[!is.na(input_vis$wTO), ]
  input_vis = plyr::arrange(input_vis, input_vis$Node.1, input_vis$Node.2)
  nodes <- data.frame(id = sort(unique(c(as.character(input_vis$Node.1), 
                                         as.character(input_vis$Node.2)))))
  g = igraph::graph_from_data_frame(input_vis, directed = F)
  DEGREE = as.data.frame(igraph::degree(g))
  igraph::E(g)$weight = abs(input_vis$wTO)
  names(DEGREE) = "degree"
  DEGREE$id = row.names(DEGREE)
  nodes = suppressMessages(plyr::join(nodes, DEGREE))
  shape.df = data.frame(shape)
  colnames(shape.df) = c("shape", "id")
  nodes=left_join(nodes, shape.df, by= "id")
  nodes$shape[is.na(nodes$shape) ] = "dot"
  nodes$value = (nodes$degree - min(nodes$degree))/(max(nodes$degree) - 
                                                      min(nodes$degree))
  nodes$value = nodes$value * 2 + 1
  nodes$size = nodes$value
  if (MakeGroups == FALSE) {
    group = 1
  }
  if (MakeGroups == "manual") {
    group = as.integer(as.factor(nodes$shape))
  }
  if (MakeGroups == "infomap") {
    group = igraph::cluster_infomap(g)$membership
  }
  else if (MakeGroups == "walktrap") {
    group = igraph::cluster_walktrap(g)$membership
  }
  else if (MakeGroups == "leading_eigen") {
    group = igraph::cluster_leading_eigen(g)$membership
  }
  else if (MakeGroups == "louvain") {
    group = igraph::cluster_louvain(g)$membership
  }
  else if (MakeGroups == "label_prop") {
    group = igraph::cluster_label_prop(g)$membership
  }
  else if (MakeGroups == "fast_greedy") {
    group = igraph::cluster_fast_greedy(g)$membership
  }
  else if (MakeGroups == "optimal") {
    group = igraph::cluster_optimal(g)$membership
  }
  else if (MakeGroups == "spinglass") {
    group = igraph::cluster_spinglass(g)$membership
  }
  else if (MakeGroups == "edge.betweenness") {
    group = igraph::edge.betweenness.community(g)$membership
  }
  nodes = plyr::join(nodes, data.frame(id = igraph::V(g)$name, 
                                       group = group))
  nodes$label = nodes$id
  nodes$title = paste0("<p> Node ID: ", nodes$id, "<br>Degree: ", 
                       nodes$degree, "</p>")
  edges <- data.frame(from = input_vis$Node.1, to = input_vis$Node.2)
  wto = abs(input_vis$wTO)
  edges$width = 0.5 + 5 * abs((wto - min(wto))/(max(wto) - 
                                                  min(wto)))
  edges$color = ifelse(input_vis$wTO > 0, "violetred", "springgreen")
  edges$title = paste0("<p> wTO: ", round(input_vis$wTO, 2), 
                       "</p>")
  ledges <- data.frame(color = c("violetred", "springgreen"), 
                       label = c("+ wTO", "- wTO"), arrows = c("", ""))
  network <- visNetwork::visNetwork(nodes, edges) %>% visNetwork::visInteraction(navigationButtons = TRUE) %>% 
    visNetwork::visEdges(smooth = smooth.edges) %>% visNetwork::visOptions(highlightNearest = list(enabled = TRUE, 
                                                                                                   degree = 1, hover = T), nodesIdSelection = list(enabled = TRUE, 
                                                                                                                                                   style = "width: 200px; height: 26px;\n   background: #f8f8f8;\n   color: darkblue;\n   border:none;\n   outline:none;"), 
                                                                           manipulation = F) %>% visNetwork::visPhysics(enabled = F) %>% 
    visNetwork::visExport(type = "pdf", name = "networkpdf", 
                          float = "left", label = "Save pdf", background = "transparent", 
                          style = "") %>% visNetwork::visExport(type = "png", 
                                                                name = "networkpng", float = "right", label = "Save png", 
                                                                background = "transparent", style = "")
  if (Cluster == T) {
    network <- network %>% visNetwork::visClusteringByGroup(groups = unique((nodes$group)))
  }
  if (legend == T) {
    network <- network %>% visNetwork::visLegend(width = 0.3, 
                                                 position = "right", main = "Group", addEdges = ledges, 
                                                 ncol = 2)
  }
  if (!is.null(layout)) {
    network <- network %>% visNetwork::visIgraphLayout(layout = layout)
  }
  if (manipulation == T) {
    network <- network %>% visNetwork::visOptions(manipulation = TRUE)
  }
  if (is.null(path)) {
    network
  }
  else if (!is.null(path)) {
    visNetwork::visSave(network, file = path)
    message(path)
  }
  nodesout = data.frame(id = nodes$id, group = nodes$group, 
                        degree = nodes$degree)
  return(list(Nodes = nodesout, network = network))
}


