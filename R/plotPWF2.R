#' @title Plot the probability Weighting Function
#' @description
#' @usage
#' plotPWF2(pwf, binsize = "auto", pwf_col = 3, pwf_lwd = 2,
#' xlab = "Biased Data in <binsize> gene bins.", ylab = "Proportion DE", ...)
#' @param pwf probability weigth function
#' @param binsize the number of gene in each bin(gene set)
#' @param pwf_col the color for the fitted line of probability weigth function
#' @param pwf_lwd the font for the fitted line of probability weigth function
#' @param xlab the label for x-axis
#' @param ylab the label for y-axis
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' png("~/GOSJ/Figure/pwfGeneGL.tiff")
#' plotPWF2(Gene.based.DE.feature.based.DE$pwfGeneGL)
#' dev.off()
#'
#' png("~/GOSJ/Figure/pwfGeneFeature.tiff")
#' plotPWF2(Gene.based.DE.feature.based.DE$pwfGeneFeature)
#' dev.off()
#'
#' png("~/GOSJ/Figure/pwfFeatureGL.tiff")
#' plotPWF2(Gene.based.DE.feature.based.DE$pwfFeatureGL)
#' dev.off()
#'
#' png("~/GOSJ/Figure/pwFeatureFeature.tiff")
#' plotPWF2(Gene.based.DE.feature.based.DE$pwfFeatureFeature)
#' dev.off()
#'
#'
plotPWF2<-function (pwf, binsize = "auto", pwf_col = 3, pwf_lwd = 2, xlab = "Biased Data in <binsize> gene bins.",
                    ylab = "Proportion DE", ...)
{
  w = !is.na(pwf$bias.data)
#  print(w)
  o = order(pwf$bias.data[w])
#  print(o)

  rang = max(pwf$pwf, na.rm = TRUE) - min(pwf$pwf, na.rm = TRUE)
  if (rang == 0 & binsize == "auto")
    binsize = 1000
  if (binsize == "auto") {
    binsize = max(1, min(100, floor(sum(w) * 0.08)))
    resid = rang
    oldwarn = options()$warn
    options(warn = -1)
    while (binsize <= floor(sum(w) * 0.1) & resid/rang >
           0.001) {
      binsize = binsize + 100
      splitter = ceiling(1:length(pwf$DEgenes[w][o])/binsize)
      de = sapply(split(pwf$DEgenes[w][o], splitter), mean)
      binlen = sapply(split(as.numeric(pwf$bias.data[w][o]),
                            splitter), mean)
      resid = sum((de - approx(pwf$bias.data[w][o], pwf$pwf[w][o],
                               binlen)$y)^2)/length(binlen)
    }
    options(warn = oldwarn)
  }
  else {
    splitter = ceiling(1:length(pwf$DEgenes[w][o])/binsize)
#   print(splitter)
    de = sapply(split(pwf$DEgenes[w][o], splitter), mean)
#    print(de)
    binlen = sapply(split(as.numeric(pwf$bias.data[w][o]),
                          splitter), median)
#    print(binlen)
  }
  xlab = gsub("<binsize>", as.character(binsize), xlab)
  if ("xlab" %in% names(list(...))) {
    if ("ylab" %in% names(list(...))) {
      plot(binlen, de, ...)
    }
    else {
      plot(binlen, de, ylab = ylab, ...)
    }
  }
  else if ("ylab" %in% names(list(...))) {
    plot(binlen, de, xlab = xlab, ...)
  }
  else {
    plot(binlen, de, xlab = xlab, ylab = ylab, ...)
  }
  lines(pwf$bias.data[w][o], pwf$pwf[w][o], col = pwf_col,
        lwd = pwf_lwd)

  return(de)

}
