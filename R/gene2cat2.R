#' Read gmt file, and convert this gmt file to a list with the name of each element being gene name
#' and each element being the pathways that this gene corresponds to
#'
#' @param gmt_input_file
#'
#' @return
#' @export
#'
#' @examples
#' gene.2.cat.hallmark<-gene2cat2("/media/H_driver/2015/Nimer_Cheng/h.all.v5.1.symbols.gmt")
#'
gene2cat2 <- function(gmt_input_file) {

  re<-GSA.read.gmt(gmt_input_file)
  gene.name<-unique(do.call(c,re$genesets))

  #gene2cat(re,"JUNB")
  #gene2cat(gene.name[1])
  gene.2.cat<-sapply(gene.name,gene2cat,re)
  #gene2cat[[]]<-re$geneset.names[res2]
  names(gene.2.cat)<-gene.name
  gene.2.cat
}
