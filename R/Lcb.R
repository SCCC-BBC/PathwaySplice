#' LrCheckBias
#'
#' Logistic regression to check bias
#'
#' @param genewise.pvalue: whether you use smallest p value or geneWisePadj
#'
#' @param sig.threshold: threshold to define differential gene list
#'
#' @param type: whether you are interested in exon or splicing junction
#'
#' @param jscs_genewise_object: an object returned from ReformatData and MakeGeneWiseTable_JunctionSeq
#'
#' @return return results from logistic regression
#' @export Boxplot of number of splicing junctions 
#'
#' @examples
#'
#' data(mds)
#' re<-Lcb(mds)
#'
Lcb <- function(jscs_genewise_object,genewise.pvalue="geneWisePadj",sig.threshold=0.05,type=c("exon","splicing")) {

  mydata<-jscs_genewise_object

  n.gene<-dim(mydata)[1]

  DE.out<-rep(0,n.gene)

  if(genewise.pvalue=="geneWisePadj"){
  de.index<-which(mydata$geneWisePadj<sig.threshold)
  }

  DE.out[de.index]<-1

  mydata.2<-cbind(mydata,DE.out)

  par(mfrow=c(1,1))

  mydata.3<-mydata.2

  mydata.3[which(mydata.3$DE.out==1),]$DE.out<-"Differential genes"

  mydata.3[which(mydata.3$DE.out==0),]$DE.out<-"Non-differential genes"

  if(type=="splicing"){
  mylogit.2 <- glm(DE.out ~ as.numeric(numKnown), data = mydata.2, family = "binomial")
  re<-summary(mylogit.2)

  pvalue<-re$coefficients[2,4]
  }else{
    mylogit.2 <- glm(DE.out ~ as.numeric(numExons), data = mydata.2, family = "binomial")
    re<-summary(mylogit.2)
    pvalue<-re$coefficients[2,4]
  }

  boxplot(unlist(mydata.3[,c(10,18)]$numExons)~unlist(mydata.3[,c(10,18)]$DE.out),ylab = "Number of exons",col="lightgray")
  text(x= 2, y= 600, labels= c("",paste0("p value from logistic regression:\n\n",pvalue)), col=c(NA,"red"))

  re<-mydata.2

  return(re)

}
