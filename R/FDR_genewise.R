#' Title
#'
#' @return
#' @export
#'
#' @examples
FDR_genewise <- function() {
  cat("tab[notZero]\n")
  print(tab[notZero])
  cat(tab[notZero],file="mm.txt")

  mm=read.table("mm.txt")
  mm=as.numeric(mm)

  cat("which(notZero)\n")
  print(which(notZero))
  cat(which(notZero),file="nn.txt")
  cat(sort(which(notZero)),file="nn2.txt")

  nn=read.table("nn.txt")
  nn=as.numeric(nn)

  cat("theta\n")
  print(theta)
  cat(theta,file="theta.txt")

  theta=read.table("theta.txt")
  theta=as.numeric(theta)

  exon.all=read.table("exon_all.txt")
  exon.all=as.numeric(exon.all)


  cat(pGene,file="pGene.txt",sep="\n",append=TRUE)

  stopifnot(all(is.finite(pGene)))

  ## Determine the thetas to be used
  theta = unique(sort(pGene))

  cat(theta,file="theta.txt",sep="\n",append=TRUE)

  pGene=read.table("pGene.txt")
  pGene=as.numeric(pGene[,1])
  theta = unique(sort(pGene))

  bins   = cut(pGene, breaks=c(-Inf, as.vector(theta)), right = TRUE, include.lowest = TRUE)

  head(bins)

  counts = tabulate(bins, nbins = nlevels(bins))

  denom  = cumsum(counts)

  #theta=read.table("theta.txt")
  #exon.all=as.numeric(ex)

  cut(c(2,3,5,1,6), breaks=c(-Inf, as.vector(c(1,3,2))), right = TRUE, include.lowest = TRUE)



  length(which(exon.all==6))

  cat(length(mm),"\t",length(nn),"\t",length(theta),"\n")



  mm2=mm[1:2]
  nn2=nn[1:2]

  mm3=c(2,3,4,5)
  nn3=c(4,5,2,3)
  theta=c(0.1,0.2,0.3,0.4,0.005,0.3)

  # 2*(1-(1-theta)^4)
  #
  # 3*(1-(1-theta)^5)
  #
  # 4*(1-(1-theta)^2)
  #
  # 5*(1-(1-theta)^3)

  numerator.test    = mapply(function(m, n) {

    #cat("m\n")
    #print(m)

    #cat("n\n")
    #print(n)

    #cat("theta","\n")
    #cat(length(theta),"\n")
    #cat(length(unique(theta)),"\n")
    #print(theta)

    re<-m * (1 - (1-theta)^n)

    #cat("re\n")
    #print(re)
    #cat(re,file="re.txt")

    re

  },
  m = mm3,
  n = nn3)
}
