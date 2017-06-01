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

# Convert the results of differential usage analysis into gene based resutls
res <- makeGeneTable(featureBasedData)
```

+ Apply logistic regression model to identify bias factor
```{r eval=TRUE, r eval=TRUE}
gene.based.table <- makeGeneTable(featureBasedData)
res <- lrTestBias(gene.based.table,loc.x=2,loc.y=150,y.lim=200,boxplot.width=0.3)
#loc.x and loc.y indicates location of p-value
```

+ Perform pathwaysplice using Canonical Pathways
```{r eval=TRUE}

dir.name <- system.file('extdata', package='PathwaySplice')

canonical.pathway.file <- '10.cp.gmt.txt'

res <- gmtGene2Cat(dir.name,canonical.pathway.file,
                   'local',genomeID='hg19')

gene.based.table <- makeGeneTable(featureBasedData)
 
res1 <- runPathwaySplice(gene.based.table,genome='hg19',
                         id='ensGene',gene2cat=res,
                         method='Wallenius')

res2 <- runPathwaySplice(gene.based.table,genome='hg19',
                         id='ensGene',gene2cat=res,
                         method='Hypergeometric')

#If you are interested in other gene sets such as Transcription Factor Targets(TFT) and hallmark gene sets from http://software.broadinstitute.org/gsea/msigdb/collections.jsp, download these gmt files, then perform analysis as the above.

```

+ Build up network based on the overlap between gene sets and visualize this network

```{r eval=TRUE}
gene.based.table <- makeGeneTable(featureBasedData)
res <- runPathwaySplice(gene.based.table,genome='hg19',
                        id='ensGene',test.cats=c('GO:BP'),
                        go.size.cut=c(5,30),
                        method='Wallenius')
output.file.dir <- "~/TestNew"
enmap <- enrichmentMap(res,n=3,similarity.threshold=0,
                       output.file.dir = output.file.dir,
                      label.vertex.by.index = TRUE)
```