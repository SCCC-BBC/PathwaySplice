#' CombineCpTftGoGeneID
#'
#' @return
#' @export
#'
#' @examples
#'
#' gene.id.all<-CombineCpTftGoGeneID("/media/H_driver/Annotation/MsigDB/c2.cp.Mouse.v5.1.symbols.gmt",
#' "/media/H_driver/Annotation/MsigDB/c3.tft.Mouse.v5.1.symbols.gmt","/media/H_driver/Annotation/GO_gene_mouse/gene_GO.mgi",
#' "/media/H_driver/Annotation/mm10/genes_table_02052016.csv")
#'
#'
CombineCpTftGoGeneID <- function(cp_input,tft_input,Go_input,gene_anno){

  gene.2.cat.cp.mouse<-Gmt2GeneCat(cp_input,gene_anno)

  gene.2.cat.tft.mouse<-Gmt2GeneCat(tft_input,gene_anno)

  gene.id.from.go<-ExtractGeneIDFromGo(Go_input,gene_anno)

  gene.id.all<-unique(c(names(gene.2.cat.cp.mouse),names(gene.2.cat.tft.mouse),gene.id.from.go))

  return(gene.id.all)

}
