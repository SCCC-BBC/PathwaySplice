#Use JunctionSeq to get analysis results
#' Get analysis results using JunctionSeq
#'
#' @param dir.name
#' @param file.sample
#' @param file.count
#' @param file.gff
#'
#' @return
#' @export
#'
#' @examples
#'
#' load("/Volumes/Bioinformatics\$/2015/Nimer_Cheng/1_29_2016.RData")
#'
#' dir.name="/Volumes/Bioinformatics$/2015/Nimer_Cheng/"
#' file.sample="decoder.bySample.txt"
#' file.count="_junction_seq_new_gtf_7/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
#' file.gff="Mus_musculus.GRCm38.83.JunctionSeq.flat.gff"
#' Re<-GetResultsFromJunctionSeq(dir.name,file.sample,file.count,file.gff)
#' head(fData(Re))
#' save(Re,file="Re_Run_test_GOSJ.RData")
#'
#' dir.name.PJ="/media/H_driver/PJ/"
#' file.sample.PJ="decoder.bySample.rtf"
#' file.count.PJ="/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
#'
#' file.gff=paste0("/media/H_driver/2015/Nimer_Cheng/",file.gff)
#' Re.PJ<-GetResultsFromJunctionSeq(dir.name.PJ,file.sample.PJ,file.count.PJ,file.gff)
#'
#' save(Re.PJ,file="/media/H_driver/PJ/PJ_jscs.RData")
#'
#'
#'buildAllPlots(jscs=jscs,outfile.prefix="./plots_based_on_DE_splice_site_gene1/",
#'gene.list=gene.based.de.splice.site[1],use.plotting.device="png",plot.gene.level.expression=TRUE,sequencing.type="single-end");
#'
#'
#'
#'
#'
GetResultsFromJunctionSeq<-function(dir.name,file.sample,file.count,file.gff){
#Get sample file
path.file.sample<-paste0(dir.name,file.sample)
decoder.bySample<-read.table(path.file.sample,header=T,stringsAsFactors = F)

print(decoder.bySample)

#Get count file
path.file.count<-paste0(dir.name,decoder.bySample$sample.ID,file.count)
countFiles<-paste0(path.file.count)

print(countFiles)

#Get annotation file
#path.file.gff<-paste0(dir.name,file.gff)
path.file.gff<-file.gff

print(path.file.gff)

#Analysis
jscs.2<-runJunctionSeqAnalyses(sample.files= countFiles,sample.names= decoder.bySample$sample.ID,condition= decoder.bySample$group.ID,
                             flat.gff.file=path.file.gff,nCores=1,verbose=TRUE,debug.mode=TRUE,use.multigene.aggregates = TRUE)

return(jscs.2)
}
