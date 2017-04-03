#' getresultsfromjunctionseq
#' 
#' This function is used to get analysis results from using JunctionSeq
#'
#' @param dir.name Path name for sample information file 
#' @param sample.file Sample information file
#' @param count.file Count file
#' @param gff.file Annotation file
#' @param method.dispFinal Determine the method 
#' used to get a 'final' dispersion estimate.
#'
#' @return The analysis result from JunctionSeq R package
#' 
#' @export
#'
#' @examples
#' 
#' dir.name <- system.file('extdata', package='PathwaySplice')
#' sample.file <- 'Sample_info.txt'
#' count.file <- 'Counts.10.genes.txt'
#' gff.file <- 'flat.chr22.10.genes.gff'
#' res <- getresultsfromjunctionseq(dir.name, sample.file, 
#' count.file,gff.file, method.dispFinal = "max")
#' 
#' 

getresultsfromjunctionseq <- function(dir.name, sample.file, 
  count.file, gff.file, method.dispFinal = c("shrink", "max", 
    "fitted", "noShare")) {
  
  # set up method for calculating dispFinal
  method.dispFinal <- match.arg(method.dispFinal)
  
  # Get sample file
  dir.name <- reformatpath(dir.name)
  
  path.sample.file <- file.path(dir.name, sample.file)
  decoder.bySample <- read.table(path.sample.file, header = TRUE, 
    stringsAsFactors = FALSE)
  
  # Get count file
  path.count.file <- file.path(dir.name, decoder.bySample$sample.ID, 
    count.file)
  
  # Get annotation file
  path.gff.file <- file.path(dir.name, "GTF_Files", gff.file)
  
  # Analysis using exonsOnly,and adjust Gender
  jscs <- runJunctionSeqAnalyses(sample.files = path.count.file, 
    sample.names = decoder.bySample$sample.ID, condition = decoder.bySample$group.ID, 
    flat.gff.file = path.gff.file, analysis.type = "exonsOnly", 
    nCores = 1, use.covars = decoder.bySample[, "Gender", 
      drop = FALSE], test.formula0 = ~sample + countbin + 
      Gender:countbin, test.formula1 = ~sample + countbin + 
      Gender:countbin + condition:countbin, effect.formula = ~condition + 
      Gender + countbin + Gender:countbin + condition:countbin, 
    geneLevel.formula = ~Gender + condition, verbose = TRUE, 
    debug.mode = TRUE, use.multigene.aggregates = TRUE, method.dispFinal = method.dispFinal)
  
  return(jscs)
}