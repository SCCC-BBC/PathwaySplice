#' Title
#'
#' @param jscs_object
#'
#' @return
#' @export
#'
#' @examples
#' re.PJ.gene.based.testable.reformat<-ReformatData(re.PJ.gene.based)
#' UseLogistic2CKBias(re.PJ.gene.based.testable.reformat)
#'
UseLogistic2CKBias <- function(jscs_object) {
  #mydata <- read.csv(input_file)

  mydata<-jscs_object
  #mydata<-mydata[-which(as.numeric(mydata$numKnown)==548),]

  ## view the first few rows of the data
  head(mydata)

  print(colnames(mydata))

  n.gene<-dim(mydata)[1]

  DE.out<-rep(0,n.gene)

  de.index<-which(mydata$geneWisePadj<0.05)

  DE.out[de.index]<-1

  mydata.2<-cbind(mydata,DE.out)

  print(head(mydata.2))

  hist(as.numeric(mydata.2$geneWisePadj))

  print(colnames(mydata.2))

  mylogit.2 <- glm(DE.out ~ as.numeric(numKnown), data = mydata.2, family = "binomial")

  logi.hist.plot(as.numeric(mydata.2$numKnown),mydata.2$DE.out,boxp=TRUE,type="hist",col="gray",xlabel="Number of splicing junctions"
                 ,counts=T)

  print(summary(mylogit.2))

}
