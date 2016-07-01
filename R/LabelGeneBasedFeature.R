#' LabelGeneBasedFeature
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
#'
#'LabelGeneBasedFeature(Re.PJ,cutoff_FC,cutoff_p,Output_file)
#'
#'
LabelGeneBasedFeature<-function(Re.PJ,cutoff_FC,cutoff_p,Output_file){

  feature.based.Re.PJ<-fData(Re.PJ)

  re2<-feature.based.Re.PJ[which(feature.based.Re.PJ[,20]>cutoff_FC&feature.based.Re.PJ[,11]<cutoff_p),]

  DE.gene<-unique(re2$geneID)

  DE_or_not<-rep(0,dim(feature.based.Re.PJ)[1])

  re3<-cbind(feature.based.Re.PJ,DE_or_not)

  re3[which(re3$geneID %in% DE.gene),]$DE_or_not<-1

  dataset2<- re3

  dataset2[sapply(dataset2, is.list)] <-
    sapply(dataset2[sapply(dataset2, is.list)],
           function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )

  write.table(dataset2,file=Output_file,row.names = FALSE,quote=FALSE,sep="\t")


  return(re3)

}
