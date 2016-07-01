#' ExtractGeneIDFromGo
#'
#' @return
#' @export
#'
#' @examples
#' gene.id.from.go<-ExtractGeneIDFromGo("/media/H_driver/Annotation/GO_gene_mouse/gene_GO.mgi","/media/H_driver/Annotation/mm10/genes_table_02052016.csv")
#'
#'
ExtractGeneIDFromGo<- function(gene_go_file,gene_anno_file) {

  GO.gene.mouse <- read.table(gene_go_file,header=F)

  colnames(GO.gene.mouse)=c("gene_id","GO")

  gene.ID.conversion<-read.csv(gene_anno_file)

  names.gene.gmt.2<-match(GO.gene.mouse$gene_id,gene.ID.conversion$gene_id)

  gene.ID.conversion.2<-gene.ID.conversion[names.gene.gmt.2,]

  #gene.2.cat.gmt.2<-gene.2.cat.gmt

  ensembl_gene_id.from.go<-unique(unlist(split(gene.ID.conversion.2[,3],";")))

  ensembl_gene_id.from.go.2<-ensembl_gene_id.from.go[-which(is.na(ensembl_gene_id.from.go))]

  ensembl_gene_id.from.go.2

}
