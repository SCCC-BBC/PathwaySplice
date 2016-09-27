InstallOtherPackage <- function() {

  #install

  #installed.packages()
  install.packages("RCurl")
  install.packages("roxygen2")
  source("http://bioconductor.org/biocLite.R")


  biocLite("org.Mm.eg.db")
  biocLite("DEXSeq")
  biocLite("biomaRt")
  biocLite("JunctionSeq")
  biocLite("goseq")
  biocLite("org.Hs.eg.db")
  biocLite("geneLenDataBase")

  install.packages("GSA")
  install.packages("VennDiagram")
  install.packages("openssl")
  biocLite("cummeRbund")
  install.packages("popbio")
  biocLite("affycoretools")
  biocLite("ReportingTools")
  biocLite("hgu95av2.db")
  install.packages("devtools")
  install.packages("rgl")
  install.packages("qpcR")
  biocLite("ggbio")

  #loading
  library(biomaRt)
  library(roxygen2)
  library(JunctionSeq)
  library(goseq)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(GO.db)
  library(geneLenDataBase)
  library(GSA)
  library(VennDiagram)
  library(cummeRbund)
  library(popbio)
  library(affycoretools)
  library(ReportingTools)
  library("hgu95av2.db")
  library(devtools)
  library(rgl)
  library(qpcR)

  # gene.model<-read.table("/media/H_driver/Annotation/mm10/genes_table_02052016.csv",header = TRUE, sep = ",", as.is=TRUE)
  # save(gene.model,file="~/GOSJ/data/gene.model.RData")

  # gene.model.hg38<-read.table("/media/H_driver/Annotation/hg38/genes_table_02092016.csv",header = TRUE, sep = ",", as.is=TRUE)
  # save(gene.model.hg38,file="~/GOSJ/data/gene.model.hg38.RData")
}

