[![Travis-CI Build Status](https://travis-ci.org/SCCC-BBC/PathwaySplice.svg?branch=master)](https://travis-ci.org/SCCC-BBC/PathwaySplice)
[![codecov](https://codecov.io/github/SCCC-BBC/PathwaySplice/coverage.svg?branch=master)](https://codecov.io/github/SCCC-BBC/PathwaySplice)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/PathwaySplice)]

# PathwaySplice
An R package for adjusting bias in pathway analysis using differential exon and splicing junction usage based results

# To Install

```{r eval=TRUE}
#In R console
library(devtools)
install_github("SCCC-BBC/PathwaySplice")

#If you use command line in pegasus terminal
R -e 'library(devtools);install_github("SCCC-BBC/PathwaySplice")'
```

# Use PathwaySplice

+ Run DEXSeq or JunctionSeq to get differential exon and/or splicing junction analysis resutls 

```{r eval=FALSE}

# test create jscs object
library(PathwaySplice)
file.sample="Sample_info.txt"
dir.name=dirname(system.file("extdata","Sample_info.txt", package = "PathwaySplice"))
dir.name=paste0(dir.name,"/")
file.gff="flat.chr22.gff"
file.count="/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
Re.example<-GetResultsFromJunctionSeq(dir.name,file.sample,file.count,file.gff)

# Convert the results of differential usage analysis into gene based resutls
library(Biobase)
re.example.gene.based<-makeGeneWiseTable(Re.example,
gene.list=unique(as.character(fData(Re.example)$geneID)))
```
+ Apply logistic regression model to identify bias factor
```{r eval=TRUE}
data(mds11)
mds33<-mds.11.sample[which(as.numeric(mds.11.sample$numExons)<=50),]
re<-LRtestBias(mds33, loc.x=2, loc.y=70,y_lim=80,boxplot_width=0.3) #loc.x and loc.y indicates location of p-value
```

+ Perform pathwaysplice in one step
```{r eval=TRUE}

data(mds11)
data(hg19)

# This function performes gene set analysis that adjusts for different number of gene subfeatures 
# (i.e. exons, spicing junctions, or gene length indicated by *adjust*) associated with each gene. 
# The significant genes are those with perGeneQvalue (from DEXSeq) < *threshold*. 

Example.Go.adjusted.by.exon<-Run_pathwaysplice(mds.11.sample,adjust='exon_SJ',sub_feature='E',
threshold=0.05,genomeID='hg19',geneID='ensGene',gene_model=hg19,method='Wallenius')

set.seed(100)
Example.Go.adjusted.by.exon.by.sampling<-Run_pathwaysplice(mds.11.sample,adjust='exon_SJ',sub_feature='E',
threshold=0.05,genomeID='hg19',geneID='ensGene',gene_model=hg19,method='Sampling')

Example.Go.unadjusted<-Run_pathwaysplice(mds.11.sample,adjust='exon_SJ',sub_feature='E',
threshold=0.05,genomeID='hg19',geneID='ensGene',gene_model=hg19,method='Hypergeometric')

```

+ If you are interested in other gene sets such as Canonical Pathways(CP),Transcription Factor Targets(TFT) and hallmark gene sets from http://software.broadinstitute.org/gsea/msigdb/collections.jsp, download these .gmt files, then perform the following steps:
```{r eval=TRUE}

cp.gmt.file=system.file("extdata","c2.cp.v5.2.symbols.gmt.txt", package = "PathwaySplice")
data(hg38)
gene.2.cat.cp.hg<-Gmt2GeneCat(cp.gmt.file,'local',gene_anno=hg38)

Example.cp.adjusted.by.exon<-Run_pathwaysplice(mds.11.sample,ad='exon_SJ',sub_feature='E',
0.05,genomeID='hg19',geneID='ensGene',gene2cat=gene.2.cat.cp.hg,gene_model=hg19,method='Wallenius')

set.seed(100)
Example.cp.adjusted.by.exon.by.sampling<-Run_pathwaysplice(mds.11.sample,ad='exon_SJ',sub_feature='E',
0.05,genomeID='hg19',geneID='ensGene',gene2cat=gene.2.cat.cp.hg,gene_model=hg19,method='Sampling')

Example.cp.unadjusted<-Run_pathwaysplice(mds.11.sample,ad='exon_SJ',sub_feature='E',
0.05,genomeID='hg19',geneID='ensGene',gene2cat=gene.2.cat.cp.hg,gene_model=hg19,method='Hypergeometric')

```

+ Build up network based on the overlap between gene sets and visualize this network

```{r eval=TRUE}
#GO
# SimilarityThreshold: Gene sets with Jaccard coefficient > similarity threshold are connected by a line,
# with thickness proportional to Jaccard coefficient    
re.w.adjusted.by.Wallenius<-enrichmentMap(Example.Go.adjusted.by.exon,n=5,SimilarityThreshold=0)
re.w.adjusted.by.sampling<-enrichmentMap(Example.Go.adjusted.by.exon.by.sampling,n=5,SimilarityThreshold=0)
re.w.unadjusted<-enrichmentMap(Example.Go.unadjusted,n=5,SimilarityThreshold=0)

#Canonical Pathways
re.w.adjusted.by.Wallenius.cp<-enrichmentMap(Example.cp.adjusted.by.exon,n=5,SimilarityThreshold=0)
re.w.adjusted.by.sampling.cp<-enrichmentMap(Example.cp.adjusted.by.exon.by.sampling,n=5,SimilarityThreshold=0)
re.w.unadjusted.cp<-enrichmentMap(Example.cp.unadjusted,n=5,SimilarityThreshold=0)
```