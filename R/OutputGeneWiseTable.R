#' Title
#'
#' @param re.PJ.gene.based
#' @param output_file
#'
#' @return
#' @export
#'
#' @examples
OutputGeneWiseTable <- function(re.PJ.gene.based,output_file) {

   no.PJ.testable.index<-which(as.character(pData(re.PJ.gene.based)$mostSigID)=="character(0)")
   re.PJ.gene.based.testable<-pData(re.PJ.gene.based)[-no.PJ.testable.index,]
   dim(re.PJ.gene.based.testable)

   cor(as.numeric(re.PJ.gene.based.testable$numKnown),
   as.numeric(re.PJ.gene.based.testable$mostSigPadjust))

  #head(re.PJ.gene.based.testable)

  table.gene.based.all.3<-re.PJ.gene.based.testable
  x <- vapply(table.gene.based.all.3$mostSigID, length, 1L) ## How many items per list element
  table.gene.based.all.3<- table.gene.based.all.3[rep(rownames(table.gene.based.all.3), x), ] ## Expand the data frame
  table.gene.based.all.3$mostSigID <- unlist(table.gene.based.all.3$mostSigID, use.names = FALSE)  ## Replace w
  #dim(table.gene.based.all.3)

  #mart = useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org", dataset="mmusculus_gene_ensembl")
  #scs.ensembl.gene.is.2.gene.symbol<-getBM(attributes=c("mgi_symbol"), filters="ensembl_gene_id", values=table.gene.based.all.3$geneID, mart=mart)
  #jscs.ensembl.gene.is.2.gene.symbol.plus.ensembl.gene.id<-cbind(scs.ensembl.gene.is.2.gene.symbol,gene.based.de.splice.site)
  write.csv(table.gene.based.all.3,file = output_file)

}
