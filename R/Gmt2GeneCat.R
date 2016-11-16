#' Gmt2GeneCat
#'
#' Read a gmt file, and return a list with its name 
#' being a gene id based on gene_anno_file, and each element of this list
#' being the pathways that this gene corresponds to
#'
#' @param gmt_input_file:input file
#' @param gene_anno_file:annotation file
#' 
#'
#' @return
#' @export
#'
#' @examples
#'
#' #gene.2.cat.cp.hg<-Gmt2GeneCat("/media/H_driver/Annotation/MsigDB/c2.cp.v5.1.symbols.gmt",
#' #"/media/H_driver/Annotation/hg38/genes_table_02092016.csv")
#'
#' #gene.2.cat.tft.hg<-Gmt2GeneCat("/media/H_driver/Annotation/MsigDB/c3.tft.v5.1.symbols.gmt",
#' #"/media/H_driver/Annotation/hg38/genes_table_02092016.csv")
#'
#' gene.2.cat.hallmark.hg<-Gmt2GeneCat("/media/H_driver/Annotation/hg38/h.all.v5.1.symbols-1.gmt",
#' "/media/H_driver/Annotation/hg38/genes_table_02092016.csv")
#'
Gmt2GeneCat <- function(gmt_input_file,gene_anno_file) {

  gene.2.cat.gmt<-gene2cat2(gmt_input_file)

  names.gene.gmt<-as.data.frame(names(gene.2.cat.gmt))
  colnames(names.gene.gmt)<-"gene_id"

  dir.name=dirname(gene_anno_file)
  dir.name=reformatPath(dir.name)
  file.name=basename(gene_anno_file)
  
  #dir(paste0(dir.name,file.name))
  
  gene_anno_file=paste0(dir.name,file.name)
  
  gene.ID.conversion<-read.csv(gene_anno_file)
  names.gene.gmt.2<-match(names.gene.gmt$gene_id,gene.ID.conversion$gene_id)
  gene.ID.conversion.2<-gene.ID.conversion[names.gene.gmt.2,]
  gene.2.cat.gmt.2<-gene.2.cat.gmt

  names(gene.2.cat.gmt.2)<-gene.ID.conversion.2[,3]

  gene.2.cat.gmt.2
}

# examples
#
# gene.2.cat.hallmark<-gene2cat2("/media/H_driver/2015/Nimer_Cheng/h.all.v5.1.symbols.gmt")
#
gene2cat2 <- function(gmt_input_file) {
  
  re<-GSA.read.gmt.2(gmt_input_file)
  gene.name<-unique(do.call(c,re$genesets))
  
  gene.2.cat<-sapply(gene.name,gene2cat,re)
  names(gene.2.cat)<-gene.name
  gene.2.cat
  
}


gene2cat <- function(gene_name,re) {
  z<-re$genesets
  res <- lapply(z, function(ch) grep(gene_name, ch))
  res2<-sapply(res, function(x) length(x) > 0)
  gene2cat<-list(re$geneset.names[res2])
  gene2cat
}


# examples
#
# re.gsa<-GSA.read.gmt.2("/media/H_driver/Annotation/MsigDB/c2.cp.Mouse.v5.1.symbols.gmt")
#

GSA.read.gmt.2<-function (filename)
{
  
  dir.name=dirname(filename)
  dir.name=reformatPath(dir.name)
  file.name=basename(filename)
  
  #dir(paste0(dir.name,file.name))
  
  filename=paste0(dir.name,file.name)
  
  
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