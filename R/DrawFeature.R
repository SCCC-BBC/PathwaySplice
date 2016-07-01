
#' DrawFeature
#'
#' @param Re.PJ
#' @param out.dir
#' @param gene_id
#'
#' @return
#' @export
#'
#' @examples
#'
#' DrawFeature(Re.PJ,output.dir.name,"ENSMUSG00000026563+ENSMUSG00000040596+ENSMUSG00000089853+ENSMUSG00000103400")
#'
#'
DrawFeature <- function(Re.PJ,out.dir,gene_id) {
buildAllPlots(Re.PJ,outfile.prefix=paste0(out.dir,gene_id,"/"),gene.list=gene_id,use.plotting.device="png",plot.gene.level.expression=TRUE)
}
