library(testthat)
test_that("test lrTestBias",
          {
            res2 <- makeGeneTable(featureBasedData)
            res3 <- lrTestBias(res2,loc.x=2,loc.y=150,y.lim=200,boxplot.width=0.3)
          })