#' OutputCatBasedPwf
#'
#' Reformat genewised results from JunctionSeq, and use this one to get enriched GO terms
#'
#' @param Re.pwf.exon.sj
#' @param gene_model
#' @param gene_2_cat
#' @param Output_file
#'
#' @return
#' @export
#'
#' @examples
#'
#' RE.cp<-OutputCatBasedPwf(Re.pwf.exon.sj,gene_model=gene.model,gene_2_cat=gene.2.cat.cp.mouse,Output_file="Cp.xls")
#'
#' RE.cft<-OutputCatBasedPwf(Re.pwf.exon.sj,gene_model=gene.model,gene_2_cat=gene.2.cat.tft.mouse,Output_file="C_tft.xls")
#'
#'
#' RE.cp.gene<-OutputCatBasedPwf(Gene.non.coding.GO$pwf,gene_model=gene.model,gene_2_cat=gene.2.cat.cp.mouse,
#' Output_file="/media/H_driver/PJ/geneGL_rm_non_coding_Cp.xls")
#'
#' RE.cft.gene<-OutputCatBasedPwf(Gene.non.coding.GO$pwf,gene_model=gene.model,gene_2_cat=gene.2.cat.tft.mouse,
#' Output_file="/media/H_driver/PJ/geneGL_rm_non_coding_tft.xls")
#'
#'

OutputCatBasedPwf<-function(Re.pwf.exon.sj,gene_model,gene_2_cat,Output_file){

  Re<-goseq2(Re.pwf.exon.sj,"mm10","ensGene",gene.model=gene_model,gene2cat=gene_2_cat,use_genes_without_cat=TRUE)

  #select GO term(10<=numInCat<=300 and BP only)

  index.select<-which(Re$numInCat>=10&Re$numInCat<=300)

  Re.select<-Re[index.select,]

  over_represented_pvalue_adjusted<-p.adjust(Re.select$over_represented_pvalue,method="BH")

  Re.select.with.adjP<-cbind(Re.select[,c(1,2)],over_represented_pvalue_adjusted,Re.select[,-c(1,2,3)])


  dataset2<- Re.select.with.adjP

  Re2<-list(EnrichedGO=Re,SelectedEnrichedGO=dataset2)

  dataset2[sapply(dataset2, is.list)] <-
    sapply(dataset2[sapply(dataset2, is.list)],
           function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )

  write.table(dataset2,file=Output_file,row.names = FALSE,quote=FALSE,sep="\t")

  return(Re2)

}
