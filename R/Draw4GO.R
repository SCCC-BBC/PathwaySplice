#' Draw4GO
#'
#' @param RE.exon.sj.all.gene
#'
#' @return
#' @export
#'
#' @examples
#' Draw4GO(RE.exon.sj.all.gene)
#'
Draw4GO <- function(RE.exon.sj.all.gene.enriched) {
  x=RE.exon.sj.all.gene.enriched[1:10,1:7]
  x$term <- factor(x$term,levels = x$term[order(x$over_represented_pvalue_adjusted,decreasing = TRUE)])
  Function<-x$term
  negative_log10p=-log10(x$over_represented_pvalue_adjusted)
  ggplot(x, aes(x=Function, y=negative_log10p,fill=factor(x$category)))+geom_bar(stat="identity")+geom_hline(yintercept = -log10(0.05))+coord_flip()
}
