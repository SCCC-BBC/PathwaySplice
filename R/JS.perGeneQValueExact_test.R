##--------------------------------------------------
## Exact computation - see methods part of the paper
##---------------------------------------------------
#' Title
#'
#' @param pGene
#' @param theta
#' @param geneSplit
#'
#' @return
#' @export
#'
#' @examples
#' source("https://bioconductor.org/biocLite.R")
#' biocLite("Biobase")
#' library(Biobase)
#'
JS.perGeneQValueExact.test = function(pGene, theta, geneSplit) {
  stopifnot(length(pGene)==length(geneSplit))

  ## Compute the numerator \sum_{i=1}^M 1-(1-theta)^{n_i}
  ## Below we first identify the summands which are the same
  ## (because they have the same n_i), then do the sum via the
  ## mapply
  numExons     = listLen(geneSplit)



  cat("numExons\n")

  cat("How many unique\n")

  cat(numExons,file="exon_all.txt","\n")

  cat(length(unique(numExons)),"\n")

  cat(unique(numExons),file="exon.txt","\n")
  cat(sort(unique(numExons)),file="exon2.txt")

  print(numExons)

  tab          = tabulate(numExons)

  cat("tab\n")
  print(tab)

  notZero      = (tab>0)

  cat("tab[notZero]\n")
  print(tab[notZero])
  cat(tab[notZero],file="mm.txt")

  cat("which(notZero)\n")
  print(which(notZero))
  cat(which(notZero),file="nn.txt")
  cat(sort(which(notZero)),file="nn2.txt")


  cat("theta\n")
  print(theta)
  cat(theta,file="theta.txt")

  numerator    = mapply(function(m, n) {

    cat("m\n")
    print(m)

    cat("n\n")
    print(n)

    cat("theta","\n")
    cat(length(theta),"\n")
    cat(length(unique(theta)),"\n")
    print(theta)

    re<-m * (1 - (1-theta)^n)

    cat("re\n")
    print(re)
    cat(re,file="re.txt")

    re

    },
                        m = tab[notZero],
                        n = which(notZero))

  print(numerator)

  cat(dim(numerator))

  numerator    = rowSums(numerator)

  ## Compute the denominator: for each value of theta, the number
  ## of genes with pGene <= theta[i].
  ## Note that in cut(..., right=TRUE), the intervals are
  ## right-closed (left open) intervals.

  bins   = cut(pGene, breaks=c(-Inf, as.vector(theta)), right = TRUE, include.lowest = TRUE)

  counts = tabulate(bins, nbins = nlevels(bins))
  denom  = cumsum(counts)
  stopifnot(denom[length(denom)]==length(pGene))

  return(numerator/denom)

}
