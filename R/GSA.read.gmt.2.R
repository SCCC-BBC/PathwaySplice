#' GSA.read.gmt.2
#' Descrption: In addition to GO terms, the GOSJ package can analyze other types of gene sets defined by users, for example, the gene sets in MSigDB database.The GSA.read.gmt.2 function takes input the gene sets in .gmt format, i.e. the first column is the pathway or gene set name, and the rest columns include gene names of the genes within the pathway. The gene names are specified in gene symbols. 
#' @param gmt filename
#'
#' @return
#' @export
#'
#' @examples
#'
#' re.gsa<-GSA.read.gmt.2("/media/H_driver/Annotation/MsigDB/c2.cp.Mouse.v5.1.symbols.gmt")
#'
#' re.ggsa.2<-GSA.read.gmt("/media/H_driver/Annotation/MsigDB/c2.cp.Mouse.v5.1.symbols.gmt")

GSA.read.gmt.2<-function (filename)
{
  a = scan(filename, what = list("", ""), sep = "\t", quote = NULL,
           fill = T, flush = T, multi.line = F)
  geneset.names = a[1][[1]]
  geneset.descriptions = a[2][[1]]
  dd = scan(filename, what = "", sep = "\t", quote = NULL)
  nn = length(geneset.names)
  n = length(dd)
  ox = rep(NA, nn)
  ii = 1
  for (i in 1:nn) {
    #cat(i)
    while ((dd[ii] != geneset.names[i]) | (dd[ii + 1] !=
                                           geneset.descriptions[i])) {
      ii = ii + 1
    }
    ox[i] = ii
    ii = ii + 1
  }
  genesets = vector("list", nn)
  for (i in 1:(nn - 1)) {
    #cat(i, fill = T)
    i1 = ox[i] + 2
    i2 = ox[i + 1] - 1
    geneset.descriptions[i] = dd[ox[i] + 1]
    genesets[[i]] = dd[i1:i2]
  }
  geneset.descriptions[nn] = dd[ox[nn] + 1]
  genesets[[nn]] = dd[(ox[nn] + 2):n]
  out = list(genesets = genesets, geneset.names = geneset.names,
             geneset.descriptions = geneset.descriptions)
  class(out) = "GSA.genesets"
  return(out)
}
