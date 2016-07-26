#' OutputGOBasedDEfromFeatures3
#'
#' Reformat genewised results from JunctionSeq, and use this one to get enriched GO terms
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
#' RE.exon.sj.all.gene<-OutputGOBasedDEfromFeatures3(re.PJ.gene.based,gene.model,"GO_exon_sj_use_all_gene.xls")
#'
#'
#'
OutputGOBasedDEfromFeatures3<-function(re.PJ.gene.based,DE_type,gene.model,Output_file_dir){

  #re2<-ReformatDataUseTable(re.PJ.gene.based)

  re2<-re.PJ.gene.based

  Re.Go.adjusted.by.exon.SJ<-GotermAnalysisUseFeatureDefineDE(re2,ad="exon_SJ",sub_feature=NULL,DE_define=DE_type,gene_model=gene.model,Output_file_dir)

  #head(Re.Go.adjusted.by.exon.SJ[[1]])

  #select GO term(10<=numInCat<=300 and BP only)

  index.select<-which(Re.Go.adjusted.by.exon.SJ[[1]]$numInCat>=10&Re.Go.adjusted.by.exon.SJ[[1]]$numInCat<=300&Re.Go.adjusted.by.exon.SJ[[1]]$ontology=="BP")

  Re.Go.adjusted.by.exon.SJ.select<-Re.Go.adjusted.by.exon.SJ[[1]][index.select,]

  over_represented_pvalue_adjusted<-p.adjust(Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue,method="BH")

  Re.Go.adjusted.by.exon.SJ.select.with.adjP<-cbind(Re.Go.adjusted.by.exon.SJ.select[,c(1,2)],over_represented_pvalue_adjusted,Re.Go.adjusted.by.exon.SJ.select[,-c(1,2,3)])


  dataset2<- Re.Go.adjusted.by.exon.SJ.select.with.adjP

  dataset2[sapply(dataset2, is.list)] <-
    sapply(dataset2[sapply(dataset2, is.list)],
           function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )

  write.table(dataset2,file=paste0(Output_file_dir,"/",DE_type,".xls"),row.names = FALSE,quote=FALSE,sep="\t")

  write.table(Re.Go.adjusted.by.exon.SJ[[2]],file=paste0(Output_file_dir,"/",DE_type,"pwf.xls"),row.names = T,quote=FALSE,sep="\t")

  return(Re.Go.adjusted.by.exon.SJ)

}
