#' calculate the p-value for GO terms
#'
#' @param num_de_incat
#' @param num_incat
#' @param num_genes
#' @param num_de
#' @param weight
#'
#' @return
#' @export
#'
#' @examples
#'
#'28 	 186 	 1852 	 125 	 1.03043
#'
#'use the gene pooled from all GO terms
#'Get_Over_under_represented_pvalue(28,186,1852,125,1.03043)
#'
#'use the genes that have matched GO terms
#'
#'Get_Over_under_represented_pvalue(28,186,1840,125,1.026767)
#'
Get_Over_under_represented_pvalue<-function(num_de_incat,num_incat,num_genes,num_de,weight){
c(dWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight)
  +pWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight,lower.tail=FALSE),
  pWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight))
}
