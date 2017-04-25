library(testthat)
test_that("test lrTestBias",
          {
            res <- lrTestBias(tiny.data,loc.x=2,loc.y=150,y.lim=200,boxplot.width=0.3)
          })