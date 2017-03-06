#' gmtGene2Cat
#'
#' Read a gene set file in GMT format, and return a list with its name
#' being a gene id, and each element of this list
#' being the pathways that this gene corresponds to
#'
#' @param dir.name Directory for the gene sets in GMT format that is located in 
#' @param pathway.file Input file for the gene sets in GMT format
#' @param file.type Indicates the gene set file in GMT format is in local or url
#' @param gene_anno_file Gene annotation file supplied as a file 
#' @param genomeID Genome ("mm10","hg19" or "hg38") to be used
#'
#' @return A list with its names being geneID, its elements being the pathways
#'
#' @export
#'
#' @examples
#' 
#' dir.name <- system.file("extdata", package="PathwaySplice")
#' canonical.pathway.file <- "c2.cp.v5.2.symbols.gmt.txt"
#' res<-gmtGene2Cat(dir.name,canonical.pathway.file,'local',genomeID="hg19")
#' 
gmtGene2Cat <- function(dir.name,pathway.file,file.type,gene_anno_file=NULL,genomeID=c("mm10","hg19","hg38")) {
  
  gmt_input_file <- file.path(dir.name,pathway.file)
  
  gene.2.cat.gmt <- gene2cat2(gmt_input_file, file.type)
  
  names.gene.gmt <- as.data.frame(names(gene.2.cat.gmt))
  colnames(names.gene.gmt) <- "gene_id"
  
  if(!is.null(gene_anno_file)){
  gene.anno.dir = dirname(gene_anno_file)
  gene.annno.dir = reformatPath(gene.anno.dir)
  file.name = basename(gene_anno_file)
  
  gene_anno_file = file.path(dir.name, file.name)
  
  gene.ID.conversion <- read.csv(gene_anno_file)
  }else{
    gene.ID.conversion=match.arg(genomeID)
  }
  
  xxx <- match2Genome(gene.ID.conversion)
  
  names.gene.gmt.2 <- match(names.gene.gmt$gene_id, xxx[,1])
  
  gene.ID.conversion.2 <- xxx[names.gene.gmt.2,]
  
  gene.2.cat.gmt.2 <- gene.2.cat.gmt
  names(gene.2.cat.gmt.2) <- gene.ID.conversion.2[,2]
  gene.2.cat.gmt.2

}

gene2cat <- function(gene_name, re) {
  z <- re$genesets
  res <- lapply(z, function(ch) grep(gene_name, ch))
  res2 <- sapply(res, function(x) length(x) > 0)
  gene2cat <- list(re$geneset.names[res2])
  gene2cat
}

gsa.read.gmt <- function(filename, type) {
  if (type != "url") {
    dir.name = dirname(filename)
    dir.name = reformatPath(dir.name)
    file.name = basename(filename)
    filename = file.path(dir.name,file.name)
  }
  
  a = scan(filename, what = list("", ""), sep = "\t", quote = NULL, 
    fill = TRUE, flush = TRUE, multi.line = FALSE)
  geneset.names = a[1][[1]]
  geneset.descriptions = a[2][[1]]
  dd = scan(filename, what = "", sep = "\t", quote = NULL)
  nn = length(geneset.names)
  n = length(dd)
  ox = rep(NA, nn)
  ii = 1
  for (i in 1:nn) {
    while ((dd[ii] != geneset.names[i]) | (dd[ii + 1] != 
      geneset.descriptions[i])) {
      ii = ii + 1
    }
    ox[i] = ii
    ii = ii + 1
  }
  genesets = vector("list", nn)
  for (i in 1:(nn - 1)) {
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

gene2cat2 <- function(gmt_input_file, file.type) {
 
  re <- gsa.read.gmt(gmt_input_file, file.type)
  gene.name <- unique(do.call(c, re$genesets))
  gene.2.cat <- sapply(gene.name, gene2cat, re)
  names(gene.2.cat) <- gene.name
  gene.2.cat
  
}