testEnrichMap <- function() {
  
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
  
  enMap2(re.cp,gene.set.type="pathway",n=4,vertex.label.font = 0.05,SimilarityThreshold=0)
  enMap2(re.cp.unadjusted,gene.set.type="pathway",n=4,vertex.label.font = 0.05,SimilarityThreshold=0)
  
}