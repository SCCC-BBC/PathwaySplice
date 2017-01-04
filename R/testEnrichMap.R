#' testPathwaySplice
#'
#' Peform one-step analysis for check bias, adjusted gene set enrichment anlysis and build network
#'
#' @param gene_based_table A gene based table converted from DEXSeq or JunctionSeq resutls
#' @param which_gene_model Gene annotation file
#'
#' @return None
#'
#' @export
#'
#' @examples
#' data(mds)
#' data(hg19)
#'
#' testPathwaySplice(mds,hg19)
#'
#'
testPathwaySplice <- function(gene_based_table, which_gene_model) {
  #Check bias using logistics regression model
  re <- LRtestBias(gene_based_table,boxplot_width = 0.3)
  
  #Analysis
  Example.Go.adjusted.by.exon <-
    Run_pathwaysplice(
      gene_based_table,
      ad = "exon_SJ",
      sub_feature = "E",
      0.05,
      genomeID = "hg19",
      geneID = "ensGene",
      gene_model = which_gene_model,
      method = "Wallenius"
    )
  Example.Go.unadjusted <-
    Run_pathwaysplice(
      gene_based_table,
      ad = "exon_SJ",
      sub_feature = "E",
      0.05,
      genomeID = "hg19",
      geneID = "ensGene",
      gene_model = which_gene_model,
      method = "Hypergeometric"
    )
  
  #Check bias
  plotPWF2(Example.Go.adjusted.by.exon[[2]])
  #gene <- names(geneList)[abs(geneList) > 2]
  
  #Construct network between gene sets
  re.w.adjusted <-
    enrichmentMap(Example.Go.adjusted.by.exon,
                  n = 5,
                  SimilarityThreshold = 0)
  re.w.unadjusted <-
    enrichmentMap(Example.Go.unadjusted,
                  n = 5,
                  SimilarityThreshold = 0)
  
}