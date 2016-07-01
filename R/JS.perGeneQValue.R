#' JS.perGeneQValue
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

#'data.perGeneQ.use.large.FC.featutre<-JS.perGeneQValue(Re.PJ.selected$pvalue,
#'Re.PJ.selected$testable,Re.PJ.selected$geneID,method = JS.perGeneQValueExact)
#'
#'length(unique(as.character(Re.PJ.selected$geneID)))
#'
#'head(Re.PJ.selected)
#'
#'
#'data.perGeneQ.use.large.FC.featutre.2<-list_to_df(data.perGeneQ.use.large.FC.featutre)
#'colnames(data.perGeneQ.use.large.FC.featutre.2)<-c("geneID","geneWisePadj_FC")
#'
#'Re.PJ.selected.2<-merge(Re.PJ.selected,data.perGeneQ.use.large.FC.featutre.2,by="geneID")
#'
#'
#'head(Re.PJ.selected.2)
#'
#'re.PJ.gene.based.selected<-merge(pData(re.PJ.gene.based),data.perGeneQ.use.large.FC.featutre.2,by="geneID")

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
