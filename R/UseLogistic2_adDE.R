#' Title
#'
#' @param input_file
#' @return
#' @export
#'
#' @examples
#'
#' UseLogistic2_adDE("/Volumes/Bioinformatics$/2015/Nimer_Cheng/
#' GeneWise_jscs3_all_with_anno_2_24_2016.csv")
#'
#'
UseLogistic2_adDE <- function(input_file) {
  mydata <- read.csv(input_file)
  ## view the first few rows of the data
  #head(mydata)

  #print(colnames(mydata))

  n.gene<-dim(mydata)[1]

  DE.out<-rep(0,n.gene)

  de.index<-which(mydata[,11]<0.05)

  DE.out[de.index]<-1

  mydata.2<-cbind(mydata,DE.out)

  #print(head(mydata.2))

  #hist(mydata.2$geneWisePadj)

  print(colnames(mydata.2))

  mylogit.sj.exon <- glm(DE.out ~ numKnown+numExons, data = mydata.2, family = "binomial")

  mylogit.sj <- glm(DE.out ~ numKnown, data = mydata.2, family = "binomial")

  mylogit.exon <- glm(DE.out ~ numExons, data = mydata.2, family = "binomial")

  #print(summary(mylogit.sj.exon))

  print(summary(mylogit.sj))

  print(summary(mylogit.exon))

}
