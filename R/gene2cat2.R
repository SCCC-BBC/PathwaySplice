#' gene2cat2
#'
#' gene2cat2 is called in Gmt2GeneCat to convert a gmt file to a list with the name of each element being gene name
#' and each element being the pathways that this gene corresponds to
#'
#' @param gmt_input_file
#'
#' @return
#' @export
#'
#' @examples
#'
#' gene.2.cat.hallmark<-gene2cat2("/media/H_driver/2015/Nimer_Cheng/h.all.v5.1.symbols.gmt")
#'
gene2cat2 <- function(gmt_input_file) {

  re<-GSA.read.gmt(gmt_input_file)
  gene.name<-unique(do.call(c,re$genesets))

  gene.2.cat<-sapply(gene.name,gene2cat,re)
  names(gene.2.cat)<-gene.name
  gene.2.cat

}
