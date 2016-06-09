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
GetMgiSymbolUsingLocalDataBase<-function(data.gene,data.set){

  gene.model<-data.set

  gene.anno<-gene.model[which(gene.model[,3]==data.gene),1]

  if(length(gene.anno)==0)
  {gene.anno="unmatched"}

  return(gene.anno)
}
