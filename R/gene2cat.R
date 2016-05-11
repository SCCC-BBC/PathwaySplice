#' Given a gene name and the object built from reading Gmt file, and find the pathways that this gene corresponds to
#'
#' @param gene_name
#' @param re
#'
#' @return
#' @export
#'
#' @examples
gene2cat <- function(gene_name,re) {

  #length(gene.name)
  #which(gene.name==gene_name)
  z<-re$genesets

  res <- lapply(z, function(ch) grep(gene_name, ch))
  res2<-sapply(res, function(x) length(x) > 0)
  #re$geneset.names[res2]

  gene2cat<-list(re$geneset.names[res2])

  #names(gene2cat)<-gene_name

  gene2cat
}
