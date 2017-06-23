[![Travis-CI Build Status](https://travis-ci.org/SCCC-BBC/PathwaySplice.svg?branch=master)](https://travis-ci.org/SCCC-BBC/PathwaySplice)
[![codecov](https://codecov.io/github/SCCC-BBC/PathwaySplice/coverage.svg?branch=master)](https://codecov.io/github/SCCC-BBC/PathwaySplice)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/PathwaySplice)]

# PathwaySplice
An R package for adjusting bias in pathway analysis using differential exon and splicing junction usage based results

# To Install

```{r eval=TRUE}
#In R console
library(devtools)

# This version of PathwaySplice can be installed if your R version is >= 3.4.0
install_github("SCCC-BBC/PathwaySplice",ref = 'development')

#In pegasus terminal 
R -e 'library(devtools);install_github("SCCC-BBC/PathwaySplice",ref="development")'

```

# Use PathwaySplice

+ Run DEXSeq or JunctionSeq to get differential exon and/or splicing junction analysis resutls 

```{r eval=TRUE}
# Load feature table to get differential usage results

library(PathwaySplice)
data(featureBasedData)

output.file.dir <- "~/OutputTestPathwaySplice"

# Convert the results of differential usage analysis into gene based resutls
gene.based.table <- makeGeneTable(featureBasedData)
```

+ Apply logistic regression model to identify bias factor
```{r eval=TRUE}
lrTestBias(gene.based.table,boxplot.width=0.3)
```

+ Perform pathwaysplice using Canonical Pathways
```{r eval=TRUE}
dir.name <- system.file('extdata', package='PathwaySplice')
#canonical.pathway.file <- '10.cp.gmt.txt'
canonical.pathway.file <- "h.all.v6.0.symbols.gmt.txt"

cpp <- gmtGene2Cat(dir.name,canonical.pathway.file,
                   'local',genomeID='hg19')

res1 <- runPathwaySplice(gene.based.table,genome='hg19',
                         id='ensGene',gene2cat=cpp,
                         method='Wallenius',go.size.limit=c(0,20))
                         
res2 <- runPathwaySplice(gene.based.table,genome='hg19',
                         id='ensGene',gene2cat=cpp,
                         method='Hypergeometric',go.size.limit=c(0,20))

output.file.name.1 <- 'In_ad_not_un_cp.xls'
output.file.name.2 <- 'In_un_not_ad_cp.xls'
compareResults(10,res1,res2,output.file.dir,
                  type.boxplot='Only3',
                  output.file.name.1,output.file.name.2)

#If you are interested in other gene sets such as Transcription Factor Targets(TFT) and hallmark gene sets from http://software.broadinstitute.org/gsea/msigdb/collections.jsp, download these gmt files, then perform analysis as the above.

dir.name.1 <- "~/Dropbox (BBSR)/Aimin_project/Research/PathwaySplice/data"
dir.name.1 <- "C:/Users/lxw391/Dropbox (BBSR)/Aimin_project/Research/PathwaySplice/data"
pathway.file <- "h.all.v6.0.symbols.gmt.txt"

cpp <- gmtGene2Cat(dir.name.1,pathway.file,'local',genomeID='hg19')

res3 <- runPathwaySplice(gene.based.table,genome='hg19',
                         id='ensGene',gene2cat=cpp,
                         method='Wallenius')
                         
res4 <- runPathwaySplice(gene.based.table,genome='hg19',
                         id='ensGene',gene2cat=cpp,
                         method='Hypergeometric')

output.file.name.3 <- 'In_ad_not_un_hp.xls'
output.file.name.4 <- 'In_un_not_ad_hp.xls'

compareResults(10,res3,res4,output.file.dir,
                  type.boxplot='Only3',
                  output.file.name.3,output.file.name.4)
                  
```

+ Build up network based on the overlap between gene sets and visualize this network

```{r eval=TRUE}
res1 <- runPathwaySplice(gene.based.table,genome='hg19',
                        id='ensGene',test.cats=c('GO:BP'),
                        go.size.limit=c(5,30),
                        method='Wallenius')
           
PathwaySplice:::writeTibble(res1,output.file.dir)
            
res2 <- runPathwaySplice(gene.based.table,genome='hg19',
                        id='ensGene',test.cats=c('GO:BP'),
                        go.size.limit=c(5,30),
                        method='Hypergeometric')
                        
enmap1 <- enrichmentMap(res1,n=6,similarity.threshold=0,
                       output.file.dir = output.file.dir,
                      label.vertex.by.index = TRUE)
                      
enmap2 <- enrichmentMap(res2,n=6,similarity.threshold=0,
                       output.file.dir = output.file.dir,
                      label.vertex.by.index = TRUE)
                      
output.file.name.1 <- 'In_ad_not_un_bp.xls'
output.file.name.2 <- 'In_un_not_ad_bp.xls'

compareResults(20,res1,res2,output.file.dir,
              type.boxplot='Only3',
              output.file.name.1,output.file.name.2)

#Label network by index of gene set, and output the network file in GML format that
#can be used as an input in Cytoscape  
enmap <- enrichmentMap(res1,n=10,similarity.threshold=0.3,label.vertex.by.index = TRUE,output.file.dir=file.path(output.file.dir,"OutEnmap"))

#Label network by description of gene set, and output the network file in GML format
#that can be used as an input in Cytoscape                       
enmap <- enrichmentMap(res1,n=10,fixed = FALSE,similarity.threshold=0.3,
                       label.vertex.by.index = FALSE,
                       output.file.dir=file.path(output.file.dir,"OutEnmap"))
```