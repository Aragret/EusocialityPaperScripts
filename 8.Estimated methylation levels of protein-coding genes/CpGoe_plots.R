setwd("/EBB-ZIVNFS/amikhail/termite_genomes/eusociality_paper/analyses/cpg_oe/prime_region")

# install.packages('LaplacesDemon')

library(LaplacesDemon)

list.files()

for (i in list.files()){
  # i = "Apac_CpGoe.txt"
  df = read.table(i, header = TRUE, sep = '\t')
  
  sp_name = substr(i, 1, 4)
  
  outfile_name = paste0(sp_name, '_CpGoe.pdf', sep = '')
  
  pdf(outfile_name)
  
  plot(density(df$CpGoe), xlim=c(0,2), lwd=2, main=sp_name, ylab="density", xlab="CpGoe")
  abline(v=Modes(df$CpGoe)$modes,lwd=2)
  
  dev.off()
}
