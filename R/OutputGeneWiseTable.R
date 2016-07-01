#' Title
#'
#' @param re.PJ.gene.based
#' @param output_file
#'
#' @return
#' @export
#'
#' @examples
#'
#' data(gene.model)
#'
#' OutputGeneWiseTable(re.PJ.gene.based,gene.model,output_file="/media/H_driver/PJ/GeneWise_Re_annotated_using_new_annotation_3.csv")
#'
#' OutputGeneWiseTable(re.PJ.gene.based,gene.model,output_file="/media/H_driver/PJ/GeneWise_Re_annotated_using_new_annotation_4.csv")
#'
#' OutputGeneWiseTable(re.PJ.gene.based,gene.model,output_file="/media/H_driver/PJ/GeneWise_41540.csv")
#'
#' OutputGeneWiseTable(Re.PJ,output_file="/media/H_driver/PJ/AllFeature_Based.csv")
#' write.table(fData(Re.PJ),row.names = FALSE,file="/media/H_driver/PJ/AllFeature_Based.csv", quote=FALSE, sep="\t")
#'
#'
OutputGeneWiseTable <- function(re.PJ.gene.based,gene.model,output_file) {

   no.PJ.testable.index<-which(as.character(pData(re.PJ.gene.based)$mostSigID)=="character(0)")
   re.PJ.gene.based.testable<-pData(re.PJ.gene.based)[-no.PJ.testable.index,]
   dim(re.PJ.gene.based.testable)

   cor(as.numeric(re.PJ.gene.based.testable$numKnown),
   as.numeric(re.PJ.gene.based.testable$mostSigPadjust))

  table.gene.based.all.3<-re.PJ.gene.based.testable
  x <- vapply(table.gene.based.all.3$mostSigID, length, 1L) ## How many items per list element
  table.gene.based.all.3<- table.gene.based.all.3[rep(rownames(table.gene.based.all.3), x), ] ## Expand the data frame
  table.gene.based.all.3$mostSigID <- unlist(table.gene.based.all.3$mostSigID, use.names = FALSE)  ## Replace w

  gene.annotation.4.geneID<-do.call(rbind,lapply(table.gene.based.all.3[,1],GenerateGeneAnno,gene.model))

  table.gene.based.all.4<-merge(gene.annotation.4.geneID,table.gene.based.all.3,by="geneID")

  write.csv(table.gene.based.all.3,file = output_file)

}

