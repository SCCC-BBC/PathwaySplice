#!/usr/bin/env Rscript
#Usage: Rscript Run_DESeq_3.R inputfile_count inputfile_sample_information outputfile_prefix

#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

#.libPaths("~/R_local_libs")

#install.packages("installr",repos="http://cran.rstudio.com/")
#library(installr)
#updateR()

#.libPaths()

#library("DESeq2")
#library("base")
#library("heatmap3")
#library("lattice")
#library("reshape")
#library("ggplot2")
#library("grid")
#library(gplots)
#library(RColorBrewer)
#library(survival)
#library(limma)
#library(edgeR)
#library(multtest)
#library(biomaRt)
#library(goseq)

#args = commandArgs(trailingOnly=T)

DE_based_on_gene_count <- function(input.file,sample.info.file,output.file) {
  # input.file = args[1]
  # sample.info.file=args[2]
  # output.file = args[3]

  cat(input.file,"\n")
  cat(sample.info.file,"\n")

  data.count<-read.table(input.file,sep="\t",header = T,row.names=1)
  data.info<-read.table(sample.info.file,sep="\t",header = F)

  print(head(data.count))
  print(data.info)

  #data.sample<-unique(trimws(as.character(data.info[,2])))
  #data.number.sample<-length(data.sample)

  #print(data.sample)
  #print(data.number.sample)


  #Function.Check.condition<-function(data.sample,data.info){


  Re<-table(data.info[,2:3])
  print(Re)

  Re1<-as.data.frame.matrix(Re)
  print(Re1)

  print(rownames(Re1))
  print(colnames(Re1))

  data.sample<-rownames(Re1)
  data.condition<-colnames(Re1)

  data.number.sample<-length(data.sample)

  sample.index<-as.list(c(seq(1,data.number.sample)))

  Function.Get.DE.FC<-function(i){

    count.4.sample.control<-data.info[which(data.info[,2] %in% data.sample[i]&trimws(data.info[,3]) %in% c(data.condition[1])),1]
    count.4.sample.treatment<-data.info[which(data.info[,2] %in% data.sample[i]&trimws(data.info[,3]) %in% c(data.condition[2])),1]

    index.control=which(as.character(gsub("X","",colnames(data.count))) %in% c(as.character(count.4.sample.control)))
    index.treatment=which(as.character(gsub("X","",colnames(data.count))) %in% c(as.character(count.4.sample.treatment)))

    data.count.reformated<-cbind(data.count[,index.control],data.count[,index.treatment])

    numControl<-length(count.4.sample.control)
    numTreatment<-length(count.4.sample.treatment)

    condition <- factor(c(rep(data.condition[1],numControl),rep(data.condition[2],numTreatment)))

    dds <- DESeqDataSetFromMatrix(data.count.reformated, DataFrame(condition), ~ condition)
    re<-results(DESeq(dds))
    re.FC<-cbind(as.data.frame(re),2^re[,2])
    colnames(re.FC)[7]="FoldChange"
    write.csv(cbind(re.FC,as.data.frame(counts(dds))),file=paste(output.file,"_",data.sample[i],"_DE.csv",sep=""))
    write.csv(as.data.frame(colData(dds)),file=paste(output.file,"_",data.sample[i],"_info.csv",sep=""))

  }

  lapply(sample.index, Function.Get.DE.FC)
}
