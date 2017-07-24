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
#output.file.dir <-"C:/Temp"

+ Perform pathwaysplice using Canonical Pathways
```{r eval=TRUE}
dir.name <- system.file('extdata', package='PathwaySplice')
#canonical.pathway.file <- '10.cp.gmt.txt'
canonical.pathway.file <- "h.all.v6.0.symbols.gmt.txt"

cpp <- gmtGene2Cat(file.path(dir.name,canonical.pathway.file),genomeID='hg19')

res1 <- runPathwaySplice(gene.based.table,genome='hg19',
                         id='ensGene',gene2cat=cpp,
                         method='Wallenius',go.size.limit=c(0,20),
                         output.file=file.path(output.file.dir,"hm_adjusted_1.csv"))
                         
res2 <- runPathwaySplice(gene.based.table,genome='hg19',
                         id='ensGene',gene2cat=cpp,
                         method='Hypergeometric',go.size.limit=c(0,20),
                         output.file=file.path(output.file.dir,"hm_unadjusted_2.csv"))

compareResults(30,res1,res2,gene.based.table,output.file.dir,
                  type.boxplot='Only3')

#If you are interested in other gene sets such as Transcription Factor Targets(TFT) and hallmark gene sets from http://software.broadinstitute.org/gsea/msigdb/collections.jsp, download these gmt files, then perform analysis as the above.

dir.name.1 <- "~/Dropbox (BBSR)/Aimin_project/Research/PathwaySplice/data"
#dir.name.1 <- "C:/Users/lxw391/Dropbox (BBSR)/Aimin_project/Research/PathwaySplice/data"
pathway.file <- "h.all.v6.0.symbols.gmt.txt"

cpp <- gmtGene2Cat(file.path(dir.name.1,pathway.file),genomeID='hg19')

res3 <- runPathwaySplice(gene.based.table,genome='hg19',
                         id='ensGene',gene2cat=cpp,
                         method='Wallenius',
                         output.file=file.path(output.file.dir,"hm_unadjusted_3.csv"))
                         
res4 <- runPathwaySplice(gene.based.table,genome='hg19',
                         id='ensGene',gene2cat=cpp,
                         method='Hypergeometric',
                         output.file=file.path(output.file.dir,"hm_unadjusted_4.csv"))

yy <- compareResults(10,res3,res4,gene.based.table,output.file.dir,
                  type.boxplot='Only3')

yy <- compareResults(10,res3,res4,gene.based.table,output.file.dir,
                  type.boxplot='All')
                  
```

+ Build up network based on the overlap between gene sets and visualize this network

```{r eval=TRUE}
res1 <- runPathwaySplice(gene.based.table,genome='hg19',
                        id='ensGene',test.cats=c('GO:BP'),
                        go.size.limit=c(5,30),
                        method='Wallenius',
                        output.file=file.path(output.file.dir,"bp_adjusted.csv"))
           
res2 <- runPathwaySplice(gene.based.table,genome='hg19',
                        id='ensGene',test.cats=c('GO:BP'),
                        go.size.limit=c(5,30),
                        method='Hypergeometric',
                        output.file=file.path(output.file.dir,"bp_unadjusted.csv"))
                        
enmap1 <- enrichmentMap(res1,n=6,similarity.threshold=0,
                       output.file.dir = output.file.dir,
                      label.node.by.index = TRUE)
                      

compareResults(20,res1,res2,output.file.dir,
              type.boxplot='Only3')

#Label network by index of gene set, and output the network file in GML format that
#can be used as an input in Cytoscape  
enmap <- enrichmentMap(res1,n=10,similarity.threshold=0.3,label.node.by.index = TRUE,output.file.dir=file.path(output.file.dir,"OutEnmap"))

#Label network by description of gene set, and output the network file in GML format
#that can be used as an input in Cytoscape                       
enmap <- enrichmentMap(res1,n=10,fixed = FALSE,similarity.threshold=0.3,
                       label.node.by.index = FALSE,
                       output.file.dir=file.path(output.file.dir,"OutEnmap"))
```

+ Some examples when using PathwaySplice on windows

```{r eval=TRUE}
# not run, demonstrates how output file can be specified
gene.based.table <- makeGeneTable(featureBasedData)

res <- runPathwaySplice(gene.based.table,genome='hg19',id='ensGene',
                          test.cats=c('GO:BP'),
                          go.size.limit=c(5,30),
                          method='Wallenius',binsize=20, 
                          output.file="C:/Users/lxw391/TEMP/test.csv")
                          
res <- runPathwaySplice(gene.based.table,genome='hg19',id='ensGene',
                          test.cats=c('GO:BP'),
                          go.size.limit=c(5,30),
                          method='Wallenius',binsize=20, 
                          output.file="C:/temp/test.csv")
```

```{r eval=TRUE}
gene.based.table <- makeGeneTable(featureBasedData)
 
res <- runPathwaySplice(gene.based.table,genome='hg19',
                          id='ensGene',test.cats=c('GO:BP'),
                          go.size.limit=c(5,30),method='Wallenius')
                          
# labeling each node by gene set name
enmap <- enrichmentMap(res,n=10,fixed = FALSE,similarity.threshold=0.3,
label.node.by.index = FALSE)
 
# labeling each node by gene set index
enmap <- enrichmentMap(res,n=10,similarity.threshold=0.3,
label.node.by.index = TRUE)
 
# not run, illustrates specification of output file directory 
enmap <- enrichmentMap(res,n=10,similarity.threshold=0.3,
label.node.by.index = TRUE, output.file.dir="C:/temp")
 
# not run, illustrates specification of output file directory
enmap <- enrichmentMap(res,n=10,similarity.threshold=0.3,
label.node.by.index = FALSE, output.file.dir="C:/temp")
```

```{r eval=TRUE}
dir.name <- system.file('extdata', package='PathwaySplice')
hallmark.pathway.file <- file.path(dir.name,"h.all.v6.0.symbols.gmt.txt")
 
hallmark <- gmtGene2Cat(hallmark.pathway.file,genomeID='hg19')
                    
gene.based.table <- makeGeneTable(featureBasedData)
 
res.adj <- runPathwaySplice(gene.based.table,genome='hg19',
                          id='ensGene',gene2cat=hallmark,  
                          go.size.limit = c(5, 200),
                          method='Wallenius', output.file=tempfile())
 
res.unadj <- runPathwaySplice(gene.based.table,genome='hg19',
                          id='ensGene',gene2cat=hallmark,go.size.limit = c(5, 200),
                          method='Hypergeometric',output.file=tempfile())
 
compareResults(20, res.adj, res.unadj, gene.based.table, type.boxplot='Only3')
 
# not run, illustrate specification of output directory
compareResults(20, res.adj, res.unadj, gene.based.table, type.boxplot='Only3',output.dir="C:/Temp")
```

+ More examples

```{r eval=TRUE}

#Exons and junctions
load("~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/peng/count_strand_based/Output_jscs/jscs.RData")

res.peng <- PathwaySplice:::makeFeatureTable(jscs)

gene.based.table.peng <- makeGeneTable(res.peng)

res.path.peng <- runPathwaySplice(gene.based.table.peng,genome='mm10',id='ensGene',test.cats=c('GO:BP'),go.size.limit=c(10,300),method='Wallenius',binsize=20,output.file=file.path("~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/peng/count_strand_based/Output_jscs","bp_adjusted.csv"))

enrichmentMap(res.path.peng,n=10,fixed = FALSE,similarity.threshold=0.3,label.node.by.index = FALSE,output.file.dir=file.path("~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/peng/count_strand_based/Output_jscs","OutEnmap"))

#Only exons
load("~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/peng/count_strand_based/Output_jscs_exonsOnly/jscs.RData")

res.peng <- PathwaySplice:::makeFeatureTable(jscs)

gene.based.table.peng <- makeGeneTable(res.peng)

res.path.peng <- runPathwaySplice(gene.based.table.peng,genome='mm10',id='ensGene',test.cats=c('GO:BP'),go.size.limit=c(10,300),method='Wallenius',binsize=20,output.file=file.path("~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/peng/count_strand_based/Output_jscsOutput_jscs_exonsOnly","bp_adjusted.csv"))

enrichmentMap(res.path.peng,n=10,fixed = FALSE,similarity.threshold=0.3,label.node.by.index = FALSE,output.file.dir=file.path("~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/peng/count_strand_based/Output_jscs_exonsOnly","OutEnmap"))

#Only junctions
load("~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/peng/count_strand_based/Output_jscs_junctionsOnly/jscs.RData")

res.peng <- PathwaySplice:::makeFeatureTable(jscs)

gene.based.table.peng <- makeGeneTable(res.peng)

res.path.peng <- runPathwaySplice(gene.based.table.peng,genome='mm10',id='ensGene',test.cats=c('GO:BP'),go.size.limit=c(10,300),method='Wallenius',binsize=20,output.file=file.path("~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/peng/count_strand_based/Output_jscs_junctionsOnly","bp_adjusted.csv"))

enrichmentMap(res.path.peng,n=10,fixed = FALSE,similarity.threshold=0.3,label.node.by.index = FALSE,output.file.dir=file.path("~/Dropbox (BBSR)/BBSR Team Folder/Aimin_Yan/peng/count_strand_based/Output_jscs_junctionsOnly","OutEnmap"))

```