library(testthat)
test_that("test lrTestBias",
          {
            res <- makeGeneTable(featureBasedData)
            res <- lrTestBias(res,loc.x=2,loc.y=150,y.lim=200,boxplot.width=0.3)
          })
test_that("test runPathwaySplice",
          {
          gene.based.table <- makeGeneTable(featureBasedData)
          res <- runPathwaySplice(gene.based.table,genome='hg19',id='ensGene',test.cats=c('GO:BP'),go.size.cut=c(5,30),method='Hypergeometric')
           })
test_that("test runPathwaySplice without size selection",
          {
            gene.based.table <- makeGeneTable(featureBasedData)
            res <- runPathwaySplice(gene.based.table,genome='hg19',id='ensGene',test.cats=c('GO:BP'),method='Hypergeometric')
          })