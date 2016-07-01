#' GeneBasedAnalysis
#'
#' This function use gene-based DE to get GO
#'
#'
#' @param Re.PJ.gene: gene-based DE results
#' @param gene_model: gene model
#' @param output_file: file name for venn
#'
#' @return
#' @export
#'
#' @examples
#'
#' Gene.non.coding.GO<-GeneBasedAnalysis(Re.PJ.gene.rm.non.coding,gene.model,"/media/H_driver/PJ/geneGL_rm_non_coding.xls")
#'
#'
GeneBasedAnalysis<-function(Re.PJ.gene,gene_model,output_file){

  re2<-Re.PJ.gene

  gene.based.DE<-rownames(re2[which(re2$padj<0.05),])

  all.gene.index<-rep(0,length(rownames(re2)))

  names(all.gene.index)<-rownames(re2)

  all.gene.index.gene.based<-all.gene.index
  all.gene.index.gene.based[which(names(all.gene.index.gene.based) %in% gene.based.DE)]=1

  Re3=nullp(all.gene.index.gene.based,"mm10","ensGene",plot.fit = FALSE)

  GO.wall.DE_interest=goseq2(Re3,"mm10","ensGene",gene.model=gene_model,use_genes_without_cat=TRUE)

  OuputGO(GO.wall.DE_interest,output_file)

  re4<-list(pwf=Re3,GO=GO.wall.DE_interest)

  return(re4)

}
