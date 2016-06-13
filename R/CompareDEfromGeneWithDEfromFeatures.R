
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


all.gene.index<-rep(0,length(re2$geneID))

names(all.gene.index)<-re2$geneID

all.gene.index.gene.based<-all.gene.index
all.gene.index.gene.based[which(names(all.gene.index.gene.based) %in% gene.based.DE)]=1

exonplussj=re2$numExons+re2$numKnown

pwf.DE_interest.gene.based.using.gene.length=nullp(all.gene.index.gene.based,"mm10","ensGene",plot.fit = FALSE)

pwf.DE_interest.gene.based.using.exonplussj=nullp(all.gene.index.gene.based,"mm10","ensGene",bias.data = exonplussj,
                                                  plot.fit = FALSE)


all.gene.index.feature.based<-all.gene.index
all.gene.index.feature.based[which(names(all.gene.index.feature.based) %in% feature.based.DE)]=1

pwf.DE_interest.feature.based.using.gene.length=nullp(all.gene.index.feature.based,"mm10","ensGene",plot.fit = FALSE)

pwf.DE_interest.feature.based.using.exonplussj=nullp(all.gene.index.feature.based,"mm10","ensGene",bias.data = exonplussj,
                                                  plot.fit = FALSE)


Re3<-list(GeneBased=gene.based.DE,FeatureBased=feature.based.DE,GenePlusFeature=re2,pwfGeneGL=pwf.DE_interest.gene.based.using.gene.length,
          pwfGeneFeature=pwf.DE_interest.gene.based.using.exonplussj,pwfFeatureGL=pwf.DE_interest.feature.based.using.gene.length,
          pwfFeatureFeature=pwf.DE_interest.feature.based.using.exonplussj)

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
