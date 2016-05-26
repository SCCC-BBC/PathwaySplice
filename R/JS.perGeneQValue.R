#' Title
#'
#' @param pvals
#' @param wTest
#' @param geneID
#' @param method
#'
#' @return
#' @export
#'
#' @examples
#'
#' gene.1<-results.with.count.data.2[which(results.with.count.data.2[,1]=="ENSMUSG00000000028+ENSMUSG00000005262"),]
#'
#' results.with.count.data.2[which(results.with.count.data.2[,1]=="ENSMUSG00000000028+ENSMUSG00000005262"),]
#'
#' sink("output.txt")
#' data.perGene<-JS.perGeneQValue(results.with.count.data.2$pvalue,results.with.count.data.2$testable,results.with.count.data.2$geneID,method = JS.perGeneQValueExact)
#' sink()
#' length(data.perGene)
#'
#'
#' data.perGene.1<-JS.perGeneQValue(gene.1$pvalue,gene.1$testable,gene.1$geneID,method = JS.perGeneQValueExact)
#' length(data.perGene.1)
#'
JS.perGeneQValue = function(pvals, wTest, geneID, method = JS.perGeneQValueExact) {

  ## use only those exons that were testable
  pvals     = pvals[wTest]
  ## 'factor' removes ununsed levels
  cat(geneID,file="geneID_0.txt",sep="\n",append=TRUE)

  cat(wTest,file="wTest.txt",sep="\n",append=TRUE)

  cat(geneID[wTest],file="geneID_wTest.txt",sep="\n",append=TRUE)

  geneID    = factor(geneID[wTest])

  cat(geneID,file="geneID.txt",sep="\n",append=TRUE)

  geneSplit = split(seq(along=geneID), geneID)

  cat(unlist(geneSplit),file="geneSplit.txt",sep="\n",append=TRUE)

  ## summarise p-values of exons for one gene: take the minimum
  pGene = sapply(geneSplit, function(i) min(pvals[i]))

  #print(pGene)
  cat(pGene,file="pGene.txt",sep="\n",append=TRUE)

  stopifnot(all(is.finite(pGene)))

  ## Determine the thetas to be used
  theta = unique(sort(pGene))

  cat(theta,file="theta.txt",sep="\n",append=TRUE)

  ## compute q-values associated with each theta
  q = method(pGene, theta, geneSplit)

  ## return a named vector of q-values per gene
  res        = rep(NA_real_, length(pGene))
  res        = q[match(pGene, theta)]
  res = pmin(1, res)
  names(res) = names(geneSplit)
  #stopifnot(!any(is.na(res)))
  return(res)
}
