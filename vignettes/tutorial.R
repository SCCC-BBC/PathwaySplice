## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(fig.width=12,fig.height=8,collapse = FALSE, comment = "#>")

## ----eval=TRUE,warning=FALSE,message=FALSE,results='hide'----------------
library(PathwaySplice)
data("featureBasedData")
gene.based.table <- makeGeneTable(featureBasedData)

## ----eval=TRUE,warning=FALSE,message=FALSE,results='hide'----------------
lrTestBias(gene.based.table,boxplot.width=0.3)

## ----eval=TRUE,warning=FALSE,message=FALSE,results='hide'----------------

# This function performes gene set analysis that adjusts for different number of gene features 
# (i.e. exons and/or spicing junctions) associated with each gene. 
# The significant genes are those with smallest p-value of gene feature  < *threshold*.

res <- runPathwaySplice(gene.based.table,genome='hg19',id='ensGene',test.cats=c('GO:BP'),go.size.limit=c(5,30),method='Wallenius')

## ----eval=TRUE,warning=FALSE,message=FALSE,results='hide'----------------
output.file.dir <- file.path(tempdir(),"OutputEnmap")
enmap <- enrichmentMap(res,n=3,output.file.dir
                       =output.file.dir,similarity.threshold=0)

