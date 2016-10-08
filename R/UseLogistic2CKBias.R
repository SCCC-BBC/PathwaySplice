#' UseLogistic2CKBias
#'
#' @param genewise.pvalue: whether you use smallest p value or geneWisePadj
#'
#' @param sig.threshold: threshold to define differential gene list
#'
#' @param type: whether you are interested in exon or splicing junction
#'
#' @param jscs_genewise_object: an object returned from ReformatData and MakeGeneWiseTable_JunctionSeq
#'
#' @return
#' @export
#'
#' @examples
#' re.PJ.gene.based.testable.reformat<-ReformatData(re.PJ.gene.based)
#' UseLogistic2CKBias(re.PJ.gene.based.testable.reformat)
#'
#' data(mds)
#' re<-UseLogistic2CKBias(mds)
#'
UseLogistic2CKBias <- function(jscs_genewise_object,genewise.pvalue="geneWisePadj",sig.threshold=0.05,type=c("exon","splicing")) {
  #mydata <- read.csv(input_file)

  mydata<-jscs_genewise_object
  #mydata<-mydata[-which(as.numeric(mydata$numKnown)==548),]

  ## view the first few rows of the data
  #head(mydata)

  #print(colnames(mydata))

  n.gene<-dim(mydata)[1]

  DE.out<-rep(0,n.gene)

  if(genewise.pvalue=="geneWisePadj"){
  de.index<-which(mydata$geneWisePadj<sig.threshold)
  }

  DE.out[de.index]<-1

  mydata.2<-cbind(mydata,DE.out)

  #print(head(mydata.2))

  par(mfrow=c(2,1))

  #hist(as.numeric(mydata.2$geneWisePadj),xlab="p value",main="Histogram")

  mydata.3<-mydata.2

  mydata.3[which(mydata.3$DE.out==1),]$DE.out<-"Differential genes"

  mydata.3[which(mydata.3$DE.out==0),]$DE.out<-"Non-differential genes"

  boxplot(unlist(mydata.3[,c(10,18)]$numExons)~unlist(mydata.3[,c(10,18)]$DE.out),ylab = "Number of exons",col="lightgray")

  #boxplot(unlist(mydata.3[,c(10,18)]$numExons)~unlist(mydata.3[,c(10,18)]$DE.out),notch=TRUE,add=TRUE,ylab = "Number of exons",col="blue")

  #print(colnames(mydata.2))

  if(type=="splicing"){
  mylogit.2 <- glm(DE.out ~ as.numeric(numKnown), data = mydata.2, family = "binomial")
  logi.hist.plot(as.numeric(mydata.2$numKnown),mydata.2$DE.out,boxp=TRUE,type="hist",col="gray",xlabel="Number of splicing junctions"
                 ,counts=T)
  a=median(as.numeric(mydata.2[which(mydata.2$DE.out==1),]$numExons))
  abline(v=a,col = "blue")
  b=median(as.numeric(mydata.2[which(mydata.2$DE.out==0),]$numExons))
  re<-summary(mylogit.2)

  pvalue<-re$coefficients[2,4]
  abline(v=b,col = "green")
  text(600,0.6,pvalue)
  }else{
    mylogit.2 <- glm(DE.out ~ as.numeric(numExons), data = mydata.2, family = "binomial")
    logi.hist.plot(as.numeric(mydata.2$numExons),mydata.2$DE.out,boxp=TRUE,type="hist",col="gray",xlabel="Number of exons"
                   ,counts=F)

    a=median(as.numeric(mydata.2[which(mydata.2$DE.out==1),]$numExons))
    abline(v=a,col = "blue")
    b=median(as.numeric(mydata.2[which(mydata.2$DE.out==0),]$numExons))
    re<-summary(mylogit.2)
    pvalue<-re$coefficients[2,4]
    abline(v=b,col = "green")
    text(600,0.6,pvalue)
  }

  re<-mydata.2

  return(re)

}
