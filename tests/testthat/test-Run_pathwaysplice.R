library(testthat)

test_that("test lrTestBias",
          {
            gene.based.table <- makeGeneTable(featureBasedData)
            lrTestBias(gene.based.table,boxplot.width=0.3)
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
test_that("test compareResults",
          {
            gene.based.table <- makeGeneTable(featureBasedData)
            
            res1 <- runPathwaySplice(gene.based.table,genome='hg19',
                                      id='ensGene',test.cats=c('GO:BP'),
                                      go.size.cut=c(10,300),method='Wallenius')
            
            res2 <- runPathwaySplice(gene.based.table,genome='hg19',
                                    id='ensGene',test.cats=c('GO:BP'),
                                    go.size.cut=c(10,300),method='Hypergeometric')
            
            output.file.dir <- "~/TestNew"
            output.file.name.1 <- 'In_ad_not_un.xls'
            output.file.name.2 <- 'In_un_not_ad.xls'

            res3 <- compareResults(25,res1,res2,output.file.dir,
                                  type.boxplot='Only3',
                                  output.file.name.1,output.file.name.2)
})

