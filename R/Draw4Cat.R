#' Draw4Cat
#'
#' @param RE.exon.sj.all.gene
#'
#' @return
#' @export
#'
#' @examples
#' Draw4Cat(RE.exon.sj.all.gene)
#'
Draw4Cat <- function(RE.exon.sj.all.gene.enriched) {
  x=RE.exon.sj.all.gene.enriched[1:10,1:5]
  x$category <- factor(x$category,levels = x$category[order(x$over_represented_pvalue_adjusted,decreasing = TRUE)])
  Function<-x$category
  negative_log10p=-log10(x$over_represented_pvalue_adjusted)
  ggplot(x, aes(x=Function, y=negative_log10p,fill=factor(x$category)))+geom_bar(stat="identity")+geom_hline(yintercept = -log10(0.05))+coord_flip()
}
