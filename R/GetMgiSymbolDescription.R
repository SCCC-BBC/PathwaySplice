
#' Title
#'
#' @param data.gene
#'
#' @return
#' @export
#'
#' @examples
#'
#'
GetMgiSymbolDescription<-function(data.gene){

  gene.anno<-as.character(paste(getBM(attributes=c("mgi_symbol","description"), filters="ensembl_gene_id",
                                      values=data.gene,
                                      mart=mart)$mgi_symbol,
                                getBM(attributes=c("mgi_symbol","description"), filters="ensembl_gene_id",
                                      values=data.gene,
                                      mart=mart)$description,sep="->"))
  if(length(gene.anno)==0)
  {gene.anno="unmatched"}

  return(gene.anno)
}

