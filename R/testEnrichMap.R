#' testPathwaySplice
#'
#' Peform one-step analysis for check bias, adjusted gene set 
#' enrichment anlysis and build network
#'
#' @param gene_based_table A gene based table converted from 
#'                        DEXSeq or JunctionSeq resutls
#' @param output.file.dir Directory for output
#' 
#' @return None
#'
#' @export
#'
#' @examples
#' 
#' dir.name <- system.file("extdata",package = "PathwaySplice")
#' output.file.dir <- file.path(dir.name,"OutputTest")
#' testpathwaysplice(tiny.data,output.file.dir = output.file.dir)
#'
#'
testpathwaysplice <- function(gene_based_table,output.file.dir) {
  #Check bias using logistics regression model
  res <- lrtestbias(gene_based_table,boxplot_width = 0.3)
  
  #Analysis
  res1 <-
    runpathwaysplice(
      gene_based_table,
      adjust = "exon_SJ",
      sub_feature = "E",
      0.05,
      genomeID = "hg19",
      geneID = "ensGene",
      method = "Wallenius"
    )
  res2 <-
    runpathwaysplice(
      gene_based_table,
      adjust = "exon_SJ",
      sub_feature = "E",
      0.05,
      genomeID = "hg19",
      geneID = "ensGene",
      method = "Hypergeometric"
    )
  
  #Construct network between gene sets
  
  
  
  res11 <-
    enrichmentmap(res1,
                  n = 5,
                  output.file.dir = output.file.dir,
                  SimilarityThreshold = 0)
  res22 <-
    enrichmentmap(res2,
                  n = 5,
                  output.file.dir = output.file.dir,
                  SimilarityThreshold = 0)
  
}