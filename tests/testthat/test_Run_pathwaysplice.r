library(PathwaySplice)
test_that("test Run_pathwaysplice",
          {
            data(mds)
            data(hg19)
            Example.Go.adjusted.by.exon<-Run_pathwaysplice(mds,ad='exon_SJ',sub_feature='E',
                                                           0.05,genomeID='hg19',geneID='ensGene',gene_model=hg19,method='Wallenius')  
          })
  


