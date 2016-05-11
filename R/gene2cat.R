#install.packages("GSA")
#library(GSA)

#genes <- c("ENSG00000124208", "ENSG00000182463", "ENSG00000124201", "ENSG00000124205", "ENSG00000124207")
#getgo(genes,'hg19','ensGene')

#data(genes)
#pwf <- nullp(genes,'hg19','ensGene')


#head(pwf)

#pvals <- goseq(pwf,'hg19','ensGene')
#head(pvals)


#gmt.input.file="/home/aiminyan/Downloads/h.all.v5.1.symbols.gmt"

#' Title
#'
#' @param gmt_input_file
#'
#' @return
#' @export
#'
#' @examples
gene2cat2 <- function(gmt_input_file) {

  re<-GSA.read.gmt(gmt_input_file)
  names(re)
  gene.name<-unique(do.call(c,re$genesets))

  gene2cat <- function(gene_name) {

    length(gene.name)
    which(gene.name==gene_name)
    z<-re$genesets

    res <- lapply(z, function(ch) grep(gene_name, ch))
    res2<-sapply(res, function(x) length(x) > 0)
    #re$geneset.names[res2]

    gene2cat<-list(re$geneset.names[res2])

    names(gene2cat)<-gene_name

    gene2cat
  }

  #gene2cat(gene.name[1])
  gene.2.cat<-sapply(gene.name,gene2cat)

  #gene2cat[[]]<-re$geneset.names[res2]
  gene.2.cat

}
#gene.2.cat.gmt<-gene2cat2(gmt.input.file)

#names(re$genesets)
#res

