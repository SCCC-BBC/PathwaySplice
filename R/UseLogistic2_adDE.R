#' Title
#'
#' @param input_file
#'
#'
#'
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'load("")
#' UseLogistic2_adDE("http://www.ats.ucla.edu/stat/data/binary.csv")
#' UseLogistic2_adDE("/media/H_driver/2015/Nimer_Cheng/GeneWise_jscs3_all_with_anno_2_24_2016.csv")
#'
#'
UseLogistic2_adDE <- function(input_file) {
  mydata <- read.csv(input_file)
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
