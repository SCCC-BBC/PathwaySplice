#' gene2cat
#'
#' Given a gene name and the object built from reading Gmt file, and find the pathways that this gene corresponds to
#'
#' @param gene_name
#' @param re
#'
#' @return
#' @export
#'
#' @examples
#'
#'
gene2cat <- function(gene_name,re) {
  z<-re$genesets
  res <- lapply(z, function(ch) grep(gene_name, ch))
  res2<-sapply(res, function(x) length(x) > 0)
  gene2cat<-list(re$geneset.names[res2])
  gene2cat
}
