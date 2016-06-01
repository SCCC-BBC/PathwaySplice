
#' Title
#'
#' @param Re
#'
#' @return
#' @export
#'
#' @examples
#'
#'
GetCCBetweenMostSigPadjustSJ <- function(Re) {
  #take some times to get this calculate doine
  re.gene.based<-makeGeneWiseTable(Re,gene.list=unique(as.character(fData(Re)$geneID)))

  no.testable.index<-which(as.character(pData(re.gene.based)$mostSigID)=="character(0)")
  re.gene.based.testable<-pData(re.gene.based)[-no.testable.index,]

  cc<-cor(as.numeric(re.gene.based.testable$numKnown),
      as.numeric(re.gene.based.testable$mostSigPadjust))
  return(cc)
}


