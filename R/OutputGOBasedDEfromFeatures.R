#' OutputGOBasedDEfromFeatures
#'
#' Use genewised results from JunctionSeq to get enriched GO terms
#'
#' @param re.PJ.gene.based
#' @param gene.model
#' @param Output_file
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' RE.exon.sj<-OutputGOBasedDEfromFeatures(re.PJ.gene.based,gene.model,"GO_exon_sj_3.xls")
#'
#' sink("test.txt")
#' RE.exon.sj<-OutputGOBasedDEfromFeatures(re.PJ.gene.based,gene.model,"GO_exon_sj_4.xls")
#' sink()
#'
OutputGOBasedDEfromFeatures<-function(re.PJ.gene.based,gene.model,Output_file){

  re<-pData(re.PJ.gene.based)

  no.re.testable.index<-which(as.character(re$mostSigID)=="character(0)")

  re2<-re[-no.re.testable.index,]

  print(names(re2))

  Re.Go.adjusted.by.exon.SJ<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(re2,ad="exon_SJ",sub_feature=NULL,0.05,gene_model=gene.model)

  head(Re.Go.adjusted.by.exon.SJ[[1]])

  #select GO term(10<=numInCat<=300 and BP only)

  index.select<-which(Re.Go.adjusted.by.exon.SJ[[1]]$numInCat>=10&Re.Go.adjusted.by.exon.SJ[[1]]$numInCat<=300&Re.Go.adjusted.by.exon.SJ[[1]]$ontology=="BP")

  Re.Go.adjusted.by.exon.SJ.select<-Re.Go.adjusted.by.exon.SJ[[1]][index.select,]

  over_represented_pvalue_adjusted<-p.adjust(Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue,method="BH")

  Re.Go.adjusted.by.exon.SJ.select.with.adjP<-cbind(Re.Go.adjusted.by.exon.SJ.select[,c(1,2)],over_represented_pvalue_adjusted,Re.Go.adjusted.by.exon.SJ.select[,-c(1,2,3)])


  dataset2<- Re.Go.adjusted.by.exon.SJ.select.with.adjP

  dataset2[sapply(dataset2, is.list)] <-
    sapply(dataset2[sapply(dataset2, is.list)],
           function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )

  write.table(dataset2,file=Output_file,row.names = FALSE,quote=FALSE,sep="\t")

  return(dataset2)

}
