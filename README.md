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
data(Feature10000BasedData)

# Convert the results of differential usage analysis into gene based resutls
gene.based.table <- makeGeneTable(Feature10000BasedData)
```

+ Apply logistic regression model to identify bias factor
```{r eval=TRUE}
res <- lrTestBias(gene.based.table,loc.x=2,loc.y=400,y.lim=500,boxplot.width=0.3)
#loc.x and loc.y indicates location of p-value
```

+ Perform pathwaysplice using Canonical Pathways
```{r eval=TRUE}
dir.name <- system.file('extdata', package='PathwaySplice')
canonical.pathway.file <- '10.cp.gmt.txt'
res <- gmtGene2Cat(dir.name,canonical.pathway.file,
                   'local',genomeID='hg19')

res1 <- runPathwaySplice(gene.based.table,genome='hg19',
                         id='ensGene',gene2cat=res,
                         method='Wallenius',go.size.cut=c(0,20))
                         
res2 <- runPathwaySplice(gene.based.table,genome='hg19',
                         id='ensGene',gene2cat=res,
                         method='Hypergeometric',go.size.cut=c(0,20))

output.dir <- "~/TestNew"
 
output.file.name.1 <- 'In_ad_not_un.xls'
output.file.name.2 <- 'In_un_not_ad.xls'
res3 <- compareResults(10,res1,res2,output.dir,output.dir,
                       type.boxplot='Only3',
                       output.file.name.1,output.file.name.2)

#If you are interested in other gene sets such as Transcription Factor Targets(TFT) and hallmark gene sets from http://software.broadinstitute.org/gsea/msigdb/collections.jsp, download these gmt files, then perform analysis as the above.

```

+ Build up network based on the overlap between gene sets and visualize this network

```{r eval=TRUE}
res1 <- runPathwaySplice(gene.based.table,genome='hg19',
                        id='ensGene',test.cats=c('GO:BP'),
                        go.size.cut=c(5,30),
                        method='Wallenius')
                        
res2 <- runPathwaySplice(gene.based.table,genome='hg19',
                        id='ensGene',test.cats=c('GO:BP'),
                        go.size.cut=c(5,30),
                        method='Hypergeometric')
                        
output.file.dir <- "~/TestNew"

enmap1 <- enrichmentMap(res1,n=6,similarity.threshold=0,
                       output.file.dir = output.file.dir,
                      label.vertex.by.index = TRUE)
                      
enmap2 <- enrichmentMap(res2,n=6,similarity.threshold=0,
                       output.file.dir = output.file.dir,
                      label.vertex.by.index = TRUE)
                      
output.file.name.1 <- 'In_ad_not_un.xls'
output.file.name.2 <- 'In_un_not_ad.xls'

res3 <- compareResults(5,res1,res2,output.dir,output.dir,
                       type.boxplot='Only3',
                       output.file.name.1,output.file.name.2)
```