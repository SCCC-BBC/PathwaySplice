library(testthat)
test_that("test lrtestbias",
          {
            data(mds11)
            mds33<-mds.11.sample[which(as.numeric(mds.11.sample$numExons)<=50),]
            re<-lrtestbias(mds33,loc.x=2,loc.y=70,y_lim=80,boxplot_width=0.3,type="splicing")
          })