#' GetResults4GeneLevel
#'
#' This function applies DESeq2 to gene level count to identify differentially expressed genes
#'
#' @param dir.name: path that count file and sample information file is in
#' @param file.sample: sample information file
#' @param file.count: count files
#'
#'
#' @return results from DESeq2
#' @export
#'
#' @examples
#' # For PJ project
#' dir.name.PJ.gene="/media/H_driver/PJ/"
#' file.sample.PJ.gene="decoder.bySample.rtf"
#' file.count.PJ.gene="/QC.geneCounts.formatted.for.DESeq.txt"
#'
#' Re.PJ.gene<-GetResults4GeneLevel(dir.name.PJ.gene,file.sample.PJ.gene,file.count.PJ.gene)

GetResults4GeneLevel<-function(dir.name,file.sample,file.count){

  suppressPackageStartupMessages(library(DESeq2))

  path.file.sample<-paste0(dir.name,file.sample)
  decoder.bySample<-read.table(path.file.sample,header=T,stringsAsFactors = F)
  #print(decoder.bySample)

  sampleCondition <- decoder.bySample$group.ID;
  sampleName <- decoder.bySample$sample.ID;


  directory <- dir.name

  sampleFiles <- paste0(decoder.bySample$sample.ID,file.count);

  sampleTable <- data.frame(sampleName = sampleName,fileName = sampleFiles,condition = sampleCondition);

  #print(sampleTable)

  dds <-  DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = directory,design = ~ condition)

  dds2 <- DESeq(dds);
  res <- results(dds2);

  return(res)

}
