#' testPathwaySplice
#'
#' @return
#' @export
#'
#' @examples
#' testPathwaySplice()
#' 
#' 
testPathwaySplice <- function() {
  
  #loading example data
  data(mds)
  
  #Check bias using logistics regression model
  re<-LRtestBias(mds)
  
  #Analysis
  Example.Go.adjusted.by.exon<-Run_pathwaysplice(mds,ad="exon_SJ",sub_feature="E",0.05,genomeID="hg19",geneID="ensGene",gene_model=hg19.gene.model,method="Wallenius")
  Example.Go.unadjusted<-Run_pathwaysplice(mds,ad="exon_SJ",sub_feature="E",0.05,genomeID="hg19",geneID="ensGene",gene_model=hg19.gene.model,method="Hypergeometric")
  
  #Check bias
  plotPWF2(Example.Go.adjusted.by.exon[[2]])
  #gene <- names(geneList)[abs(geneList) > 2]
  
  #Construct network between gene sets
  re.w.adjusted<-enrichmentMap(Example.Go.adjusted.by.exon,n=5,SimilarityThreshold=0)
  re.w.unadjusted<-enrichmentMap(Example.Go.unadjusted,n=5,SimilarityThreshold=0)
  
}
