#  GetResultsFromJunctionSeq
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
#' # For example data set
#'
#' dir.name="/media/H_driver/Aimin_project/"
#'
#' file.sample="decoder.bySample.Mut_WT_2.rtf"
#' file.gff="Homo_sapiens.GRCh38.84.processed.sorted.4.JunctionSeq.flat.gff"
#' file.count="/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
#'
#' file.gff=paste0(dir.name,"GTF_Files/",file.gff)
#'
#' Re.example<-GetResultsFromJunctionSeq(dir.name,file.sample,file.count,file.gff)
#' save(Re.example,file=paste0(getwd(),"/data/example_Mut_WT.RData"))
#'
#' file.sample.2="decoder.bySample.Mut_HC_2.rtf"
#' Re.example.Mut.HC<-GetResultsFromJunctionSeq(dir.name,file.sample.2,file.count,file.gff)
#'
#' buildAllPlots(jscs=jscs,outfile.prefix="./plots_based_on_DE_splice_site_gene1/",
#' gene.list=gene.based.de.splice.site[1],use.plotting.device="png",plot.gene.level.expression=TRUE,sequencing.type="single-end");
#'
#'
#
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
#jscs.2<-runJunctionSeqAnalyses(sample.files= countFiles,sample.names= decoder.bySample$sample.ID,condition= decoder.bySample$group.ID,
#                             flat.gff.file=path.file.gff,nCores=1,verbose=TRUE,debug.mode=TRUE,use.multigene.aggregates = TRUE)

#Analysis using exonsOnly , and adjust Gender
jscs.2<-runJunctionSeqAnalyses(sample.files= countFiles,sample.names= decoder.bySample$sample.ID,
                               condition= decoder.bySample$group.ID,
                               flat.gff.file=path.file.gff,
                               analysis.type = "exonsOnly",nCores=1,
                               use.covars = decoder.bySample[,"Gender",drop=F],
                               test.formula0  = ~ sample + countbin + Gender : countbin,
                               test.formula1  = ~ sample + countbin + Gender : countbin + condition : countbin,
                               effect.formula = ~ condition + Gender + countbin + Gender : countbin + condition : countbin,
                               geneLevel.formula = ~ Gender + condition,
                               verbose=TRUE,debug.mode=TRUE,use.multigene.aggregates = TRUE)

return(jscs.2)
}
