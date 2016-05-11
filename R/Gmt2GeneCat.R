#' Read a gmt file, return a list with the name of element being a gene id based on gene_anno_file, and each element
#' being the pathways that this gene corresponds to
#'
#' @param gmt_input_file
#' @param gene_anno_file
#'
#' @return
#' @export
#'
#' @examples
Gmt2GeneCat <- function(gmt_input_file,gene_anno_file) {

  gene.2.cat.gmt<-gene2cat2(gmt_input_file)
  #length(gene.2.cat.gmt)

  names.gene.gmt<-as.data.frame(names(gene.2.cat.gmt))
  colnames(names.gene.gmt)<-"gene_id"
  human.gene.ID.conversion<-read.csv(gene_anno_file)
  names.gene.gmt.2<-match(names.gene.gmt$gene_id,human.gene.ID.conversion$gene_id)
  human.gene.ID.conversion.2<-human.gene.ID.conversion[names.gene.gmt.2,]
  gene.2.cat.gmt.2<-gene.2.cat.gmt
  names(gene.2.cat.gmt.2)<-human.gene.ID.conversion.2[,3]
  gene.2.cat.gmt.2
}
