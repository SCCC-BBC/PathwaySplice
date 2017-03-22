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
install_github("SCCC-BBC/PathwaySplice")

#In pegasus terminal 
R -e 'library(devtools);install_github("SCCC-BBC/PathwaySplice")'

```

# Use PathwaySplice

+ Run DEXSeq or JunctionSeq to get differential exon and/or splicing junction analysis resutls 

```{r eval=FALSE}

# test create jscs object
library(PathwaySplice)
dir.name <- system.file("extdata",package = "PathwaySplice")
sample.file <- "Sample_info.txt"
count.file <- "QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
gff.file <- "flat.chr22.gff"
res <- getresultsfromjunctionseq(dir.name,sample.file,count.file,gff.file)

# Convert the results of differential usage analysis into gene based resutls
res1 <- makegenewisetable(res)
```
+ Apply logistic regression model to identify bias factor
```{r eval=TRUE, r eval=TRUE}
res <- mds.11.sample[which(as.numeric(mds.11.sample$numExons)<=50),]
res1 <- lrtestbias(res, loc.x=2, loc.y=70,y_lim=80,boxplot_width=0.3) #loc.x and loc.y indicates location of p-value
```

+ Perform pathwaysplice using Canonical Pathways
```{r eval=TRUE}

dir.name <- system.file("extdata", package="PathwaySplice")
canonical.pathway.file <- "c2.cp.v5.2.symbols.gmt.txt"

res <- gmtgene2cat(dir.name,canonical.pathway.file,'local',genomeID="hg19")

res1 <- runpathwaysplice(mds.11.sample,adjust='exon_SJ',sub_feature='E',
0.05,genomeID='hg19',geneID='ensGene',gene2cat=res,method='Wallenius')

set.seed(100)
res2 <- runpathwaysplice(mds.11.sample,adjust='exon_SJ',sub_feature='E',
0.05,genomeID='hg19',geneID='ensGene',gene2cat=res,method='Sampling')

res3 <- runpathwaysplice(mds.11.sample,adjust='exon_SJ',sub_feature='E',
0.05,genomeID='hg19',geneID='ensGene',gene2cat=res,method='Hypergeometric')

#If you are interested in other gene sets such as Transcription Factor Targets(TFT) and hallmark gene sets from http://software.broadinstitute.org/gsea/msigdb/collections.jsp, download these gmt files, then perform analysis as the above.

```

+ Build up network based on the overlap between gene sets and visualize this network

```{r eval=TRUE}

#Canonical Pathways
res11 <- enrichmentmap(res1,n=5,SimilarityThreshold=0)
res22 <- enrichmentmap(res2,n=5,SimilarityThreshold=0)
res33 <- enrichmentmap(res3,n=5,SimilarityThreshold=0)
```