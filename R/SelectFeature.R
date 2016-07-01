#' SelectFeature
#'
#' @param Re.PJ
#'
#' @return
#' @export
#'
#' @examples
#' Re.PJ.selected<-SelectFeature(Re.PJ)
#'
#'dim(Re.PJ.selected)

SelectFeature<-function(Re.PJ){

  table.based.RE.PJ<-fData(Re.PJ)

  re2<-table.based.RE.PJ[table.based.RE.PJ$testable,]

  re2<-table.based.RE.PJ[which(table.based.RE.PJ[,20]>=1),]

  return(re2)

}
