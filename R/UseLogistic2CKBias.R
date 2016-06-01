#' Title
#'
#' @param jscs_object
#'
#' @return
#' @export
#'
#' @examples
#'
#'
UseLogistic2CKBias <- function(jscs_object) {
  #mydata <- read.csv(input_file)

  mydata<-jscs_object
  ## view the first few rows of the data
  head(mydata)

  print(colnames(mydata))

  n.gene<-dim(mydata)[1]

  DE.out<-rep(0,n.gene)

  de.index<-which(mydata[,9]<0.05)

  DE.out[de.index]<-1

  mydata.2<-cbind(mydata,DE.out)

  print(head(mydata.2))

  hist(mydata.2$geneWisePadj)

  print(colnames(mydata.2))

  mylogit.2 <- glm(DE.out ~ numKnown, data = mydata.2, family = "binomial")

  print(summary(mylogit.2))

}
