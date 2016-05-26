#' Title
#'
#' @param input_file
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
#'
#' Re<-UseOtherDataSet("/media/H_driver/DataSet_SJ/GSE66793_cmp1-3-geo-juncs.tsv")
#'
#'
UseOtherDataSet<-function(input_file){
  #load(input_file)

  Re.Jun<-read.table(input_file,sep="\t",header = TRUE)



  # Re.DE<-read.csv(gene_based_input_file)
  #
  print(colnames(Re.Jun))

  cat(dim(Re.Jun)[1],"\n")

  cat(length(unique(Re.Jun[,13])),"\n")


  # print(head(Re.DE))
  #
  #
  # Re.Jun.2<-Re.Jun[,-1]
  #
  # colnames(Re.Jun.2)[1]="ensembl_gene_id"
  #
  # colnames(Re.DE)[1]="geneID"
  #
  # mart = useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org", dataset="mmusculus_gene_ensembl")
  # geneMaps<-getBM(attributes=c("ensembl_gene_id","mgi_symbol"), filters="mgi_symbol",values=Re.DE[,1],mart=mart)
  #
  # colnames(geneMaps)=c("ensembl_gene_id","geneID")
  #
  # Re.DE.2<-merge(Re.DE,geneMaps,by="geneID")
  #
  # print(head(Re.DE.2))
  #
  # Re.DE.Jun<-merge(Re.DE.2,Re.Jun.2,by="ensembl_gene_id")
  #
  #
  # print(dim(Re.DE.Jun))
  #
  # Re.DE.Jun.2<-Re.DE.Jun[-which(is.na(Re.DE.Jun[,8])),]
  #
  # n<-dim(Re.DE.Jun.2)[1]
  # ran.p<-runif(n)
  #
  #
  # Re<-list(Re.DE.Jun.with.NA=Re.DE.Jun,Re.DE.Jun.without.NA=cbind(Re.DE.Jun.2,ran.p))
  #
  # #plot()
  #
  # pairs(~padj+geneWisePadj+mostSigPadjust+numKnown+ran.p,data=Re[[2]][,c(8,22,24,26,33)],main="4 p values vs SJ")
  # #boxplot(cbind(Re[[2]][,c(8,22,24)],ran.p))
  #
  #
  # Re

}
