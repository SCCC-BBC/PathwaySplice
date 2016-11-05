
#' testPathwaySplice
#'
#' @return
#' @export
#'
#' @examples
testPathwaySplice <- function() {
  
  library(org.Hs.eg.db)
  library(PathwaySplice)
  library(goseq)
  library(gdata)
  library(GO.db)
  library(DOSE,quietly = TRUE)
  library(reshape2)
  library(igraph)
  
  data(mds)
  
  #Analysis
  Example.Go.adjusted.by.exon<-Run_pathwaysplice(mds,ad="exon_SJ",sub_feature="E",0.05,genomeID="hg19",geneID="ensGene",gene_model=hg19.gene.model,method="Wallenius")
  
  Example.Go.unadjusted<-Run_pathwaysplice(mds,ad="exon_SJ",sub_feature="E",0.05,genomeID="hg19",geneID="ensGene",gene_model=hg19.gene.model,method="Hypergeometric")
  
  #Check bias
  plotPWF2(Example.Go.adjusted.by.exon[[2]])
  #gene <- names(geneList)[abs(geneList) > 2]
  
  #Construct network between gene sets
  re.w.adjusted<-enMap2(Example.Go.adjusted.by.exon,n=4,SimilarityThreshold=0)
  re.w.unadjusted<-enMap2(Example.Go.unadjusted,n=4,SimilarityThreshold=0)
  
  re.cp<-pathwaysplice(Example.Go.adjusted.by.exon[[2]],genome,id,gene.model=gene.model,gene2cat=gene.2.cat.cp.hg,method="Wallenius")
  re.cp.unadjusted<-pathwaysplice(Example.Go.adjusted.by.exon[[2]],genome,id,gene.model=gene.model,gene2cat=gene.2.cat.cp.hg,method="Hypergeometric")
  
  re.tft<-pathwaysplice(Example.Go.adjusted.by.exon[[2]],genome,id,gene.model=gene.model,gene2cat=gene.2.cat.tft.hg,method="Wallenius")
  re.tft.unadjusted<-pathwaysplice(Example.Go.adjusted.by.exon[[2]],genome,id,gene.model=gene.model,gene2cat=gene.2.cat.tft.hg,method="Hypergeometric")
  
  re.hallmark<-pathwaysplice(Example.Go.adjusted.by.exon[[2]],genome,id,gene.model=gene.model,gene2cat=gene.2.cat.hallmark.hg,method="Wallenius")
  re.hallmark.unadjusted<-pathwaysplice(Example.Go.adjusted.by.exon[[2]],genome,id,gene.model=gene.model,gene2cat=gene.2.cat.hallmark.hg,method="Hypergeometric")
  
  re.pathway.enrichMap.adjusted<-enMap2(re.cp,gene.set.type="pathway",n=4,vertex.label.font = 0.05,SimilarityThreshold=0)
  re.pathway.enrichMap.unadjusted<-enMap2(re.cp.unadjusted,gene.set.type="pathway",n=4,vertex.label.font = 0.05,SimilarityThreshold=0)
  
  re.tft.enrichMap.adjusted<-enMap2(re.tft,gene.set.type="pathway",n=4,vertex.label.font = 0.05,SimilarityThreshold=0)
  re.tft.enrichMap.unadjusted<-enMap2(re.tft.unadjusted,gene.set.type="pathway",n=4,vertex.label.font = 0.05,SimilarityThreshold=0)
  
  re.hallmark.enrichMap.adjusted<-enMap2(re.hallmark,gene.set.type="pathway",n=4,vertex.label.font = 0.05,SimilarityThreshold=0)
  re.hallmark.enrichMap.unadjusted<-enMap2(re.hallmark.unadjusted,gene.set.type="pathway",n=4,vertex.label.font = 0.05,SimilarityThreshold=0)
  
}
