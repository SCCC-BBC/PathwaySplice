
#' CompareDEfromGeneWithDEfromFeatures
#'
#' This function compares gene-based DE with feature-based DE
#'
#'
#' @param Re.PJ.gene: gene-based DE results
#' @param re.PJ.gene.based: feature-based DE results
#' @param Output_venn_file: file name for venn
#'
#' @return
#' @export
#'
#' @examples
#'
#' Gene.based.DE.feature.based.DE<-CompareDEfromGeneWithDEfromFeatures(Re.PJ.gene,re.PJ.gene.based,"/media/H_driver/PJ/gene_feature_DE_2.tiff")
#'
CompareDEfromGeneWithDEfromFeatures<-function(Re.PJ.gene,re.PJ.gene.based,Output_venn_file){

re<-merge(Re.PJ.gene,pData(re.PJ.gene.based),by=0)

no.re.testable.index<-which(as.character(re$mostSigID)=="character(0)")

re2<-re[-no.re.testable.index,]


gene.based.DE<-re2[which(re2$padj<0.05),]$geneID

feature.based.DE<-re2[which(re2$geneWisePadj<0.05),]$geneID

Re3<-list(GeneBased=gene.based.DE,FeatureBased=feature.based.DE)

venn.plot <- venn.diagram(
  x = Re3[c(1,2)],
  filename = Output_venn_file,
  height = 3000,
  width = 3500,
  resolution = 1000,
  col = "black",
  lty = "dotted",
  lwd = 1,
  fill = c("red","blue"),
  alpha = 0.50,
  label.col = c(rep("white",3)),
  cex = 0.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("red","blue"),
  cat.cex = 0.5,
  cat.pos = 0.5,
  cat.dist = 0.05,
  cat.fontfamily = "serif"
)

return(Re3)

}

