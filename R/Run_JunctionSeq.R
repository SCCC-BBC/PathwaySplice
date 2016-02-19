#Install required packages
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("DESeq2")
biocLite("DEXSeq")
biocLite("JunctionSeq")

biocValid("JunctionSeq")
library(edgeR)
library(DESeq2)
library(DEXSeq)
library(JunctionSeq)
library(GO.db)

install.packages("fBasics")
library(fBasics)
source("http://bioconductor.org/biocLite.R")
biocLite("org.Mm.eg.db")
library(org.Mm.eg.db)
install.packages("VennDiagram")
library(VennDiagram)
require("gplots")
sessionInfo()

#Use JunctionSeq to get analysis results
GetResultsFromJunctionSeq<-function(dir.name,file.sample,file.count,file.gff){
#Get sample file 
path.file.sample<-paste0(dir.name,file.sample)
decoder.bySample<-read.table(path.file.sample,header=T,stringsAsFactors = F)

#Get count file
path.file.count<-paste0(dir.name,decoder.bySample$sample.ID,file.count)
countFiles<-paste0(path.file.count)

#Get annotation file
path.file.gff<-file.gff

#Analysis
jscs.2<-runJunctionSeqAnalyses(sample.files= countFiles,sample.names= decoder.bySample$sample.ID,condition= decoder.bySample$group.ID,
                             flat.gff.file=path.file.gff,nCores=1,verbose=TRUE,debug.mode=TRUE,use.multigene.aggregates = TRUE)

return(jscs.2)
}

#GetResultsFromJunctionSeq("/media/H_driver/2015/Nimer_Cheng/","decoder.bySample.txt","_junction_seq_new_gtf_7/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt",
#"/media/DATA/mus_musculus/Mus_musculus.GRCm38.83.JunctionSeq.flat.gff")


