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
The latest version installed 
by visiting Bioconductor website http://bioconductor.org/packages/devel/bioc/html/PathwaySplice.html

or by 
```{r eval=FALSE, message=FALSE, warning=FALSE, results='hide'}
library(devtools)
install_github("SCCC-BBC/PathwaySplice",ref = 'development')
```
After installation, the PathwaySplice package can be loaded into R using:
```{r eval=TRUE, message=FALSE, warning=FALSE, results='hide'}
library(PathwaySplice)
```
### Reference
Yan A, Ban Y, Gao Z, Chen X, Wang L. 2017. PathwaySplice: An R package for unbiased pathway analysis of alternative splicing in RNA-Seq data. _Submitted_

Anders, Simon, Alejandro Reyes, and Wolfgang Huber. 2012. “Detecting differential usage of exons from RNA-seq data.” Genome Research 22 (10): 2008–17. doi:10.1101/gr.133744.111.

Hartley, Stephen W., and James C. Mullikin. 2016. “Detection and visualization of differential splicing in RNA-Seq data with JunctionSeq.” Nucleic Acids Research, June. Oxford University Press, gkw501. doi:10.1093/nar/gkw501.

Shannon, P., Andrew Markiel, Owen Ozier, Nitin S Baliga, Jonathan T Wang, Daniel Ramage, Nada Amin, Benno Schwikowski, and Trey Ideker. 2003. “Cytoscape: A Software Environment for Integrated Models of Biomolecular Interaction Networks.” Genome Research 13 (11): 2498–2504. doi:10.1101/gr.1239303.

Young, Matthew D., Matthew J. Wakefield, Gordon K. Smyth, and Alicia Oshlack. 2010. “Gene Ontology Analysis for Rna-Seq: Accounting for Selection Bias.” Genome Biology 11 (2): R14. doi:10.1186/gb-2010-11-2-r14.

<!-- Usage: rmarkdown::render("vignettes/tutorial.Rmd", output_format="all") --> 
<!-- Usage: rmarkdown::render("vignettes/tutorial.Rmd", output_format="all",encoding="utf8")(on windows) -->
