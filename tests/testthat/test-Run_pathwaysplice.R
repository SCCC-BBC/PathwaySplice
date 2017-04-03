library(testthat)
test_that("test lrtestbias",
          {
            res <- lrtestbias(tiny.data,loc.x=2,loc.y=150,y_lim=200,boxplot_width=0.3)
          })