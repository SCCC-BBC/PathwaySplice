#' lrtestbias
#'
#' Logistic regression to check bias
#'
#' @param jscs_genewise_object Gene based table 
#' @param genewise.pvalue Whether you use smallest p value or geneWisePadj
#' @param sig.threshold Threshold to define differential gene list
#' @param type Whether you are interested in exon or splicing junction
#' @param loc.x Indicates x coordinate of p value from logistic regression 
#' @param loc.y Indicates y coordinate of p value from logistic regression
#' @param y_lim Defining the largest number of exons in y axis in boxplot 
#' @param boxplot_width Parameter for boxplot width
#'
#' @return Results from logistic regression
#' @export
#'
#' @examples
#' res <- lrtestbias(tiny.data,loc.x=2,loc.y=150,y_lim=200,boxplot_width=0.3)

lrtestbias <- function(jscs_genewise_object, genewise.pvalue = "geneWisePadj", 
    sig.threshold = 0.05, type = c("exon","splicing"),loc.x=2,loc.y=70,y_lim=80,boxplot_width) {

    mydata <- jscs_genewise_object
    
    n.gene <- dim(mydata)[1]
    
    DE.out <- rep(0, n.gene)
    
    if (genewise.pvalue == "geneWisePadj") {
        de.index <- which(mydata$geneWisePadj < sig.threshold)
    }
    
    DE.out[de.index] <- 1
    
    mydata.2 <- cbind(mydata, DE.out)
    
    par(mfrow = c(1, 1))
    
    mydata.3 <- mydata.2
    
    mydata.3[which(mydata.3$DE.out == 1), ]$DE.out <- "Significant genes"
    
    mydata.3[which(mydata.3$DE.out == 0), ]$DE.out <- "Non significant genes"
    
    type <- match.arg(type)
    
    switch (type,
            splicing = {
              cat("Use splicing junctions\n")
              
              if(var(as.numeric(unlist(mydata.2$numKnown)))!=0){
              
              mylogit.2 <- glm(DE.out ~ as.numeric(numKnown), data = mydata.2, 
                               family = "binomial")
              re <- summary(mylogit.2)
              pvalue <- re$coefficients[2, 4]
              }else{
                cat("There are no variations on the number of splicing junctions\n")
            }
            },
            {
              cat("Use exons\n")
              mylogit.2 <- glm(DE.out ~ as.numeric(numExons), data = mydata.2, 
                               family = "binomial")
              re <- summary(mylogit.2)
              pvalue <- re$coefficients[2, 4]
              pvalue <- format(pvalue,width=8,digits=4)
              
              temp <- data.frame(mydata.3[, c(10, 18)])
              
              temp$DE.out <- factor(temp$DE.out)
              
              temp$DE.out <- factor(temp$DE.out,levels=levels(temp$DE.out)[c(2,1)])
              
              boxplot(unlist(temp$numExons) ~ unlist(temp$DE.out),
                      boxwex=boxplot_width,ylab = "Number of exons",col = "lightgray",ylim=c(1,y_lim))
              
              text(x = loc.x, y = loc.y, labels = c("", paste0("P value from logistic regression:\n\n", 
                                                               pvalue)), col = c(NA, "black"))
              
            }
    )
    
    re <- mydata.3
    
    return(re)
    
}