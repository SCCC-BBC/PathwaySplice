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
#' UseLogistic2_adDE("http://www.ats.ucla.edu/stat/data/binary.csv")
#'
UseLogistic2_adDE <- function(input_file) {
  mydata <- read.csv(input_file)
  ## view the first few rows of the data
  head(mydata)
  mydata$rank <- factor(mydata$rank)
  mylogit <- glm(admit ~ gre + gpa + rank, data = mydata, family = "binomial")
  p.re<-predict(mylogit,type="response")

  re<-cbind(mydata$admit,p.re)

  print(head(re))

}
