#' CompareGOResults
#'
#' @param Gene.based.DE.feature.based.DE
#' @param gene_model
#'
#' @return
#' @export
#'
#' @examples
#'
#'
CompareGOResults <- function(Gene.based.DE.feature.based.DE, gene_model)
  {
  GO.wall.DE_interest.geneGL=goseq2(Gene.based.DE.feature.based.DE$pwfGeneGL,"mm10","ensGene",gene.model=gene_model)

  GO.wall.DE_interes.geneFT=goseq2(Gene.based.DE.feature.based.DE$pwfGeneFeature,"mm10","ensGene",gene.model=gene_model)

  GO.wall.DE_interest.FtFT=goseq2(Gene.based.DE.feature.based.DE$pwfFeatureFeature,"mm10","ensGene",gene.model=gene_model)

  GO.wall.DE_interest.FtFT=goseq2(Gene.based.DE.feature.based.DE$pwfFeatureFeature,"mm10","ensGene",gene.model=gene_model)
}
