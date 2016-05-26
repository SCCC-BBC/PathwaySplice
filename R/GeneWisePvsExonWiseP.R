#' @title Compare gene based p value with exon based p value
#'
#' @param input_data_set data set for all genes
#' @param gene_wise_p_c genewise p value column
#' @param feature_wise_p_c feature(exon or SJ) p value column
#'
#' @return
#' @export
#'
#' @examples
#'
#' load("/media/H_driver/2015/Nimer_Cheng/Data_set_two_methods.RData")
#' GeneWisePvsExonWiseP(data.table.gene.based.all.jscs4,7,9)
#'
#'
#' Re<-GeneWisePvsExonWiseP("/media/H_driver/2015/Nimer_Cheng/GeneWise_jscs3_all_with_anno_2_24_2016.csv",
#'"/media/H_driver/2015/Nimer_Cheng/DE_cheng_output_sample1_DE.csv")
#'
#'
#'
#'
#'
#'
GeneWisePvsExonWiseP<-function(subfeature_based_input_file,gene_based_input_file){
#load(input_file)

  Re.Jun<-read.csv(subfeature_based_input_file)
  Re.DE<-read.csv(gene_based_input_file)

  print(head(Re.Jun))
  print(head(Re.DE))


  Re.Jun.2<-Re.Jun[,-1]

  colnames(Re.Jun.2)[1]="ensembl_gene_id"

  colnames(Re.DE)[1]="geneID"

  mart = useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org", dataset="mmusculus_gene_ensembl")
  geneMaps<-getBM(attributes=c("ensembl_gene_id","mgi_symbol"), filters="mgi_symbol",values=Re.DE[,1],mart=mart)

  colnames(geneMaps)=c("ensembl_gene_id","geneID")

  Re.DE.2<-merge(Re.DE,geneMaps,by="geneID")

  print(head(Re.DE.2))

  Re.DE.Jun<-merge(Re.DE.2,Re.Jun.2,by="ensembl_gene_id")


  print(dim(Re.DE.Jun))

  Re.DE.Jun.2<-Re.DE.Jun[-which(is.na(Re.DE.Jun[,8])),]

  n<-dim(Re.DE.Jun.2)[1]
  ran.p<-runif(n)


  Re<-list(Re.DE.Jun.with.NA=Re.DE.Jun,Re.DE.Jun.without.NA=cbind(Re.DE.Jun.2,ran.p))

  #plot()

  pairs(~padj+geneWisePadj+mostSigPadjust+numKnown+ran.p,data=Re[[2]][,c(8,22,24,26,33)],main="4 p values vs SJ")
  #boxplot(cbind(Re[[2]][,c(8,22,24)],ran.p))


  Re


}
