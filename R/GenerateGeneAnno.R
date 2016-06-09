#' @ Generate gene annotations based on the Ensembl ID of genes
#'
#' @param TobeAnno: gene list
#'
#' @return
#' @export
#'
#' @examples
#'
#' Cheng.gene.all.anno.3<-do.call(rbind,lapply(re.PJ.gene.based.testable[,1],GenerateGeneAnno,gene.model))
#'
#'
GenerateGeneAnno<-function(TobeAnno,data.set){
  #First element
  gene.split<-unlist(strsplit(TobeAnno,"\\+"))
  #print(gene.split)

  #use bioMart data base
  #gene.anno<-GetMgiSymbolDescription(gene.split[1])

  #use local file
  gene.anno<-GetMgiSymbolUsingLocalDataBase(gene.split[1],data.set)

  if(length(gene.split)>=2){
    for(j in 2:length(gene.split)){
      #gene.anno.temp=GetMgiSymbolDescription(gene.split[j])
      gene.anno.temp<-GetMgiSymbolUsingLocalDataBase(gene.split[j],data.set)
      gene.anno<-paste(gene.anno,gene.anno.temp,sep="+")
    }
  }

  Annotatedgene<-as.data.frame(cbind(TobeAnno,gene.anno))
  colnames(Annotatedgene)=c("geneID","geneAnno")

  return(Annotatedgene)

}
