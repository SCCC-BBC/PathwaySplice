#' Gmt2GeneCat
#'
#' Read a gmt file, and return a list with the name of element being a gene id based on gene_anno_file, and each element
#' being the pathways that this gene corresponds to
#'
#' @param gmt_input_file
#' @param gene_anno_file
#' @param based_by
#'
#' @return
#' @export
#'
#' @examples
#'
#' gene.2.cat.cp.mouse<-Gmt2GeneCat("/media/H_driver/Annotation/MsigDB/c2.cp.Mouse.v5.1.symbols.gmt",
#' "/media/H_driver/Annotation/mm10/genes_table_02052016.csv")
#'
#' gene.2.cat.tft.mouse<-Gmt2GeneCat("/media/H_driver/Annotation/MsigDB/c3.tft.Mouse.v5.1.symbols.gmt",
#' "/media/H_driver/Annotation/mm10/genes_table_02052016.csv")
#'

Gmt2GeneCat <- function(gmt_input_file,gene_anno_file) {

  gene.2.cat.gmt<-gene2cat2(gmt_input_file)

  names.gene.gmt<-as.data.frame(names(gene.2.cat.gmt))
  colnames(names.gene.gmt)<-"gene_id"

  human.gene.ID.conversion<-read.csv(gene_anno_file)
  names.gene.gmt.2<-match(names.gene.gmt$gene_id,human.gene.ID.conversion$gene_id)
  human.gene.ID.conversion.2<-human.gene.ID.conversion[names.gene.gmt.2,]
  gene.2.cat.gmt.2<-gene.2.cat.gmt

  names(gene.2.cat.gmt.2)<-human.gene.ID.conversion.2[,3]

  gene.2.cat.gmt.2
}
