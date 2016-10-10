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
  biocLite("clusterProfiler")

  source("http://bioconductor.org/biocLite.R")
  biocLite("AnnotationHub")
  biocLite("Homo.sapiens")
  biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
  biocLite("BSgenome.Hsapiens.UCSC.hg19")
  biocLite("biomaRt")
  biocLite("TxDb.Athaliana.BioMart.plantsmart22")
  biocLite("OrderedList")
  
  #loading
  library(DEXSeq)
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
  library(clusterProfiler)
  library(OrderedList)

  # source("http://bioconductor.org/biocLite.R")
  # biocLite("AnnotationHub")
  # biocLite("Homo.sapiens")
  # biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
  # biocLite("BSgenome.Hsapiens.UCSC.hg19")
  # biocLite("biomaRt")
  # biocLite("TxDb.Athaliana.BioMart.plantsmart22")
  # 
  # listMarts()
  # 
  # ensembl <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")
  # ensembl
  # head(listAttributes(ensembl),20)
  # 
  # head(listDatasets(ensembl),20)
  # 
  # res <- getBM(attributes=c("hgnc_symbol","entrezgene","ensembl_gene_id","chromosome_name","start_position","end_position"), mart = ensembl)
  # 
  # hg19.gene.model=res
  # 
  # save(hg19.gene.model,file=paste0(getwd(),"/data/hg19.RData"))
  # 
  # gene.model<-read.table("/media/H_driver/Annotation/mm10/genes_table_02052016.csv",header = TRUE, sep = ",", as.is=TRUE)
  # save(gene.model,file="~/GOSJ/data/gene.model.RData")

  # gene.model.hg38<-read.table("/media/H_driver/Annotation/hg38/genes_table_02092016.csv",header = TRUE, sep = ",", as.is=TRUE)
  # save(gene.model.hg38,file="~/GOSJ/data/gene.model.hg38.RData")
}

