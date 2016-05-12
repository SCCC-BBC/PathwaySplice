#' codes for simulation and plot PWF
#' simulation code from xiaobei
#' PWF plot code from goseq
#' Get data set for simulation studies
#'
#' @param counts
#' @param drop.extreme.dispersion
#' @param drop.low.lambda
#'
#' @return
#' @export
#'
#' @examples
getDataset <- function(counts, drop.extreme.dispersion = 0.1, drop.low.lambda = TRUE) {
  ## this function generates NB parameters from real dataset ##
  ## it is low-level function of NBsim ##
  d <- DGEList(counts)
  d <- calcNormFactors(d)
  cp <- round(cpm(d,normalized.lib.sizes=TRUE),1)
  if(drop.low.lambda) d <- d[rowSums(cp>1) >= 2, ]
  d$AveLogCPM <- log2(rowMeans(cpm(d, prior.count = 1e-5)))
  d <- estimateGLMCommonDisp(d)
  d <- estimateGLMTrendedDisp(d)
  d <- estimateGLMTagwiseDisp(d)
  dispersion <- d$tagwise.dispersion
  AveLogCPM <- d$AveLogCPM
  if(is.numeric(drop.extreme.dispersion))
  {
    bad <- quantile(dispersion, 1-drop.extreme.dispersion, names = FALSE)
    ids <- dispersion <= bad
    AveLogCPM <- AveLogCPM[ids]
    dispersion <- dispersion[ids]
  }
  dataset.AveLogCPM <- AveLogCPM
  dataset.dispersion <- dispersion
  dataset.lib.size <- d$samples$lib.size
  dataset.nTags <- nrow(d)
  list(dataset.AveLogCPM = dataset.AveLogCPM, dataset.dispersion = dataset.dispersion, dataset.lib.size = dataset.lib.size, dataset.nTags = dataset.nTags)
}
