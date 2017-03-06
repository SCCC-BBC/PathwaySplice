#' GetResultsFromJunctionSeq
#' Get analysis results using JunctionSeq  #### I guess we could add a brief introduction of the results, e.g. peak results, mappting results etc
#'
#' @param dir.name path name for sample information file 
#' @param file.sample sample information file
#' @param file.count count file      #### short reads count? or other count?
#' @param file.gff annotation file
#'
#' @return return an analysis results from JunctionSeq

#' @export
#'
#' @examples
#' file.sample="Sample_info.txt"
#' dir.name=dirname(system.file("extdata","Sample_info.txt", package = "PathwaySplice"))
#' dir.name=paste0(dir.name,"/")
#' file.gff="flat.chr22.gff"
#' file.count="/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
#' \donttest{Re.example<-GetResultsFromJunctionSeq(dir.name,file.sample,file.count,file.gff)}
#' 
GetResultsFromJunctionSeq <- function(dir.name, file.sample, 
  file.count, file.gff) {

  # Get sample file
  dir.name = reformatPath(dir.name)
  
  path.file.sample <- paste0(dir.name, file.sample)
  decoder.bySample <- read.table(path.file.sample, header = TRUE, 
    stringsAsFactors = FALSE)
  
  # Get count file
  path.file.count <- paste0(dir.name, decoder.bySample$sample.ID, 
    file.count)
  countFiles <- paste0(path.file.count)
  
  # Get annotation file
  path.file.gff <- paste0(dir.name, "GTF_Files/", file.gff)
  
  # Analysis using exonsOnly , and adjust Gender
  jscs.2 <- runJunctionSeqAnalyses(sample.files = countFiles, 
    sample.names = decoder.bySample$sample.ID, condition = decoder.bySample$group.ID, 
    flat.gff.file = path.file.gff, analysis.type = "exonsOnly", 
    nCores = 1, use.covars = decoder.bySample[, "Gender", 
      drop = FALSE], test.formula0 = ~sample + countbin + 
      Gender:countbin, test.formula1 = ~sample + countbin + 
      Gender:countbin + condition:countbin, effect.formula = ~condition + 
      Gender + countbin + Gender:countbin + condition:countbin, 
    geneLevel.formula = ~Gender + condition, verbose = TRUE, 
    debug.mode = TRUE, use.multigene.aggregates = TRUE)
  
  return(jscs.2)
}
