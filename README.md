[![Travis-CI Build Status](https://travis-ci.org/SCCC-BBC/PathwaySplice.svg?branch=master)](https://travis-ci.org/SCCC-BBC/PathwaySplice)
[![codecov](https://codecov.io/github/SCCC-BBC/PathwaySplice/coverage.svg?branch=master)](https://codecov.io/github/SCCC-BBC/PathwaySplice)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/PathwaySplice)]

# PathwaySplice
Pathway analysis for alternative splicing in RNA-seq datasets that accounts for different number of gene features

## Introduction

In alternative splicing ananlysis of RNASeq data, one popular approach is to first identify gene features (e.g. exons or junctions) significantly associated with splicing using methods such as DEXSeq [Anders2012] or JunctionSeq [Hartley2016], and then perform pathway analysis based on the list of genes associated with the significant gene features. 

For DEXSeq results, we use _gene features_ to refers to non-overlapping exon counting bins [Anders2012, Figure 1], while for JunctionSeq results, _gene features_ refers to non-overlapping exon or splicing junction counting bins. 

A major challenge is that without explicit adjustment, pathways analysis would be biased toward pathways that include genes with a large number of gene features, because these genes are more likely to be selected as "significant genes" in pathway analysis.  

PathwaySplice is an R package that falicitate the folowing analysis: 

1. Performing pathway analysis that explicitly adjusts for the number of exons or junctions associated with each gene; 
2. Visualizing selection bias due to different number of exons or junctions for each gene and formally tests for presence of bias using logistic regression; 
3. Supporting gene sets based on the Gene Ontology terms, as well as more broadly defined gene sets (e.g. MSigDB) or user defined gene sets; 
4. Identifing the significant genes driving pathway significance and 
5. Organizing significant pathways with an enrichment map, where pathways with large number of overlapping genes are grouped together in a network graph.


### Installation

After installation, the PathwaySplice package can be loaded into R using:
```{r eval=TRUE, message=FALSE, warning=FALSE, results='hide'}
library(PathwaySplice)
```

The latest version can also be installed by 
```{r eval=FALSE, message=FALSE, warning=FALSE, results='hide'}
library(devtools)
install_github("SCCC-BBC/PathwaySplice",ref = 'development')
```

### Vignette

![pdf version](https://github.com/SCCC-BBC/PathwaySplice/blob/master/tutorial.pdf)


### Reference
Yan A, Ban Y, Gao Z, Chen X, Wang L (2017) PathwaySplice: An R package for unbiased pathway analysis of alternative splicing in RNA-Seq data. _Submitted_
Anders, S., Reyes, A. and Huber, W. (2012) Detecting differential usage of exons from RNA-seq data, Genome research, 22, 2008-2017.
Hartley, S.W. and Mullikin, J.C. (2016) Detection and visualization of differential splicing in RNA-Seq data with JunctionSeq, Nucleic acids research, 44, e127
<!-- Usage: rmarkdown::render("vignettes/tutorial.Rmd", output_format="all") --> 
<!-- Usage: rmarkdown::render("vignettes/tutorial.Rmd", output_format="all",encoding="utf8")(on windows) -->
