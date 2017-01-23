[![Travis-CI Build Status](https://travis-ci.org/SCCC-BBC/PathwaySplice.svg?branch=master)](https://travis-ci.org/SCCC-BBC/PathwaySplice)
[![codecov](https://codecov.io/github/SCCC-BBC/PathwaySplice/coverage.svg?branch=master)](https://codecov.io/github/SCCC-BBC/PathwaySplice)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/PathwaySplice)](https://cran.r-project.org/package=PathwaySplice)

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

# Make a tiny example data set using chr4 data
```{r eval=FALSE}

dir.name="/media/H_driver/Aimin_project/GOSJ_STAR_Bam/"
file.name=dir(dir.name,recursive = TRUE,pattern="sorted.bam")
file.name.whole<-paste0(dir.name,file.name)
file.name.selected<-file.name.whole
file.name.selected.2<-as.list(file.name.selected)
names(file.name.selected.2)=sapply(strsplit(file.name.selected,split="\\/"),"[[",6)
gtf1="/media/H_driver/Aimin_project/GTF_Files/Homo_sapiens.GRCh38.84.processed.sorted.2.chr4.only.gtf"
cmd4="java -Xmx5000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --noGzipOutput --dropChrom chr1,chr10,chr11,chr11_KI270721v1_random,chr12,chr13,chr14,chr14_GL000009v2_random,chr14_GL000194v1_random,chr14_GL000225v1_random,chr14_KI270722v1_random,chr14_KI270723v1_random,chr14_KI270726v1_random,chr15,chr15_KI270727v1_random,chr16,chr16_KI270728v1_random,chr17,chr17_GL000205v2_random,chr17_KI270729v1_random,chr18,chr19,chr1_KI270706v1_random,chr1_KI270707v1_random,chr1_KI270709v1_random,chr1_KI270710v1_random,chr1_KI270711v1_random,chr1_KI270712v1_random,chr1_KI270713v1_random,chr1_KI270714v1_random,chr2,chr20,chr21,chr22,chr22_KI270731v1_random,chr22_KI270732v1_random,chr22_KI270733v1_random,chr22_KI270734v1_random,chr22_KI270735v1_random,chr22_KI270736v1_random,chr22_KI270738v1_random,chr2_KI270716v1_random,chr3,chr3_GL000221v1_random,chr4_GL000008v2_random,chr5,chr6,chr7,chr8,chr9,chr9_KI270718v1_random,chr9_KI270719v1_random,chr9_KI270720v1_random,chrM,chrUn_GL000195v1,chrUn_GL000213v1,chrUn_GL000214v1,chrUn_GL000216v2,chrUn_GL000218v1,chrUn_GL000219v1,chrUn_GL000220v1,chrUn_GL000224v1,chrUn_KI270311v1,chrUn_KI270315v1,chrUn_KI270330v1,chrUn_KI270337v1,chrUn_KI270362v1,chrUn_KI270435v1,chrUn_KI270438v1,chrUn_KI270442v1,chrUn_KI270467v1,chrUn_KI270511v1,chrUn_KI270519v1,chrUn_KI270522v1,chrUn_KI270590v1,chrUn_KI270741v1,chrUn_KI270742v1,chrUn_KI270743v1,chrUn_KI270744v1,chrUn_KI270745v1,chrUn_KI270746v1,chrUn_KI270747v1,chrUn_KI270748v1,chrUn_KI270750v1,chrUn_KI270751v1,chrUn_KI270754v1,chrX,chrY"

re.out<-lapply(file.name.selected.2,callQoRT,output.dir="chr4_drop_other",gtf_file=gtf1,runing_cmd=cmd4)

add "Homo_sapiens.GRCh38.84.processed.sorted.4.JunctionSeq.flat.chr4.100.gff" into PathwaySplice/inst/extdata/GTF_Files/


```


# Use PathwaySplice

+ Run DEXSeq or JunctionSeq to get differential exon and/or splicing junction analysis resutls 

```{r eval=FALSE}
library(PathwaySplice)
dir.name=dirname(system.file("extdata","decoder.bySample.Mut_WT_example.txt", package = "PathwaySplice"))
dir.name=paste0(dir.name,"/")
file.sample="decoder.bySample.Mut_WT_example.txt"
file.gff="Homo_sapiens.GRCh38.84.processed.sorted.4.JunctionSeq.flat.chr4.100.gff"

file.count="/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"

Re.example<-GetResultsFromJunctionSeq(dir.name,file.sample,file.count,file.gff)
 
```

+ Convert the results of differential usage analysis into gene based resutls

```{r eval=FALSE}

all.gene.list<-unique(as.character(fData(Re.example)$geneID))
 
#make a tiny example data set
choosed.gene.list<-all.gene.list
re.example.gene.based<-makeGeneWiseTable(Re.example,
gene.list=choosed.gene.list)
tiny.data<-re.example.gene.based

```
+ Apply logistic regression model to identify bias factor
```{r eval=TRUE}
data(mds11)
mds33<-mds.11.sample[which(as.numeric(mds.11.sample$numExons)<=50),]
re<-LRtestBias(mds33,p.x=2,p.y=70,y_lim=80,boxplot_width=0.3)
```

+ Perform pathwaysplice in one step
```{r eval=TRUE}

data(mds11)
data(hg19)

Example.Go.adjusted.by.exon<-Run_pathwaysplice(mds.11.sample,ad='exon_SJ',sub_feature='E',
0.05,genomeID='hg19',geneID='ensGene',gene_model=hg19,method='Wallenius')
temp=Example.Go.adjusted.by.exon$GO.selected
mean(temp[which(temp$category %in% c("GO:0072331","GO:0072332","GO:0097193")),]$rank.value.by.over_represented_pvalue)

set.seed(100)
Example.Go.adjusted.by.exon.by.sampling<-Run_pathwaysplice(mds.11.sample,ad='exon_SJ',sub_feature='E',
0.05,genomeID='hg19',geneID='ensGene',gene_model=hg19,method='Sampling')
temp=Example.Go.adjusted.by.exon.by.sampling$GO.selected
mean(temp[which(temp$category %in% c("GO:0072331","GO:0072332","GO:0097193")),]$rank.value.by.over_represented_pvalue)

Example.Go.unadjusted<-Run_pathwaysplice(mds.11.sample,ad='exon_SJ',sub_feature='E',
0.05,genomeID='hg19',geneID='ensGene',gene_model=hg19,method='Hypergeometric')
temp=Example.Go.unadjusted$GO.selected
mean(temp[which(temp$category %in% c("GO:0072331","GO:0072332","GO:0097193")),]$rank.value.by.over_represented_pvalue)

```

+ If you are interested in other gene sets such as Canonical Pathways(CP),Transcription Factor Targets(TFT) and hallmark gene sets from http://software.broadinstitute.org/gsea/msigdb/collections.jsp, download these .gmt files, then perform the following steps:
```{r eval=TRUE}

cp.gmt.file=system.file("extdata","c2.cp.v5.2.symbols.gmt.txt", package = "PathwaySplice")
data(hg38)
gene.2.cat.hallmark.hg<-Gmt2GeneCat(cp.gmt.file,'local',gene_anno=hg38)

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
re.w.adjusted.by.Wallenius<-enrichmentMap(Example.Go.adjusted.by.exon,n=5,SimilarityThreshold=0)
re.w.adjusted.by.sampling<-enrichmentMap(Example.Go.adjusted.by.exon.by.sampling,n=5,SimilarityThreshold=0)
re.w.unadjusted<-enrichmentMap(Example.Go.unadjusted,n=5,SimilarityThreshold=0)
```
