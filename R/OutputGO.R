#' OutputGO
#'
#' @param Re.Go.adjusted.by.exon.SJ
#' @param Output_file
#'
#' @return
#' @export
#'
#' @examples
#' 
#' Output.GO<-OuputGO(Example.Go.adjusted.by.exon[[1]]$GO,Output_file=paste0(getwd(),"/example.xls"))
#' 
#' Draw4GO(Output.GO)
#'
OuputGO<- function(Re.Go.adjusted.by.exon.SJ,Output_file) {

  index.select<-which(Re.Go.adjusted.by.exon.SJ$numInCat>=10&Re.Go.adjusted.by.exon.SJ$numInCat<=300&Re.Go.adjusted.by.exon.SJ$ontology=="BP")

  Re.Go.adjusted.by.exon.SJ.select<-Re.Go.adjusted.by.exon.SJ[index.select,]

  over_represented_pvalue_adjusted<-p.adjust(Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue,method="BH")

  Re.Go.adjusted.by.exon.SJ.select.with.adjP<-cbind(Re.Go.adjusted.by.exon.SJ.select[,c(1,2)],over_represented_pvalue_adjusted,Re.Go.adjusted.by.exon.SJ.select[,-c(1,2,3)])

  dataset2<- Re.Go.adjusted.by.exon.SJ.select.with.adjP

  dataset2[sapply(dataset2, is.list)] <-
    sapply(dataset2[sapply(dataset2, is.list)],
           function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )

  #Draw4GO(dataset2)
  
  write.table(dataset2,file=Output_file,row.names = FALSE,quote=FALSE,sep="\t")
 
  return(dataset2)
}
