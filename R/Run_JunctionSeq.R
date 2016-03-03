# #Install required packages
# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
# biocLite("DESeq2")
# biocLite("DEXSeq")
# #biocLite("JunctionSeq")
#
# biocValid("JunctionSeq")
# library(edgeR)
# library(DESeq2)
# library(DEXSeq)
# library(JunctionSeq)
# library(GO.db)
#
# install.packages("fBasics")
# library(fBasics)
# biocLite("org.Mm.eg.db")
# library(org.Mm.eg.db)
# install.packages("VennDiagram")
# library(VennDiagram)
# require("gplots")
# sessionInfo()
#
# #Use JunctionSeq to get analysis results
# GetResultsFromJunctionSeq<-function(dir.name,file.sample,file.count,file.gff){
# #Get sample file
# path.file.sample<-paste0(dir.name,file.sample)
# decoder.bySample<-read.table(path.file.sample,header=T,stringsAsFactors = F)
#
# #Get count file
# path.file.count<-paste0(dir.name,decoder.bySample$sample.ID,file.count)
# countFiles<-paste0(path.file.count)
#
# #Get annotation file
# path.file.gff<-file.gff
#
# #Analysis
# jscs.2<-runJunctionSeqAnalyses(sample.files= countFiles,sample.names= decoder.bySample$sample.ID,condition= decoder.bySample$group.ID,
#                              flat.gff.file=path.file.gff,nCores=1,verbose=TRUE,debug.mode=TRUE,use.multigene.aggregates = TRUE)
#
# return(jscs.2)
# }
#
# #Run for the data from cheng
# jscs.3<-GetResultsFromJunctionSeq("/media/H_driver/2015/Nimer_Cheng/","decoder.bySample.txt","_junction_seq_new_gtf_7/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt",
#                           "/media/DATA/mus_musculus/Mus_musculus.GRCm38.83.JunctionSeq.flat.gff")
# data.jscs.3<-fData(jscs.3)
# gene.id.all.3<-unique(as.character(data.jscs.3[,1]))
# table.gene.based.all.jscs3<-makeGeneWiseTable(jscs.3,gene.id.all.3,FDR.threshold = 0.05,verbose = TRUE, debug.mode = FALSE)
# data.table.gene.based.all.jscs3<-pData(table.gene.based.all.jscs3)
#
# #extract genes that are testable
# data.table.gene.based.all.jscs4<-data.table.gene.based.all.jscs3[which(as.character(data.table.gene.based.all.jscs3[,8])!="character(0)"),]
#
# #Annotate all genes
# Cheng.gene.all.anno.new<-do.call(rbind,lapply(data.table.gene.based.all.jscs4[,1],GenerateGeneAnno))
# Cheng.gene.all.anno.2.data.table.2<-merge(Cheng.gene.all.anno.3,data.table.gene.based.all.jscs4,by="geneID",sort=FALSE)
#
# #Output results
# Cheng.gene.all.anno.2.data.table.3<- data.frame(lapply(Cheng.gene.all.anno.2.data.table.2, as.character), stringsAsFactors=FALSE)
# write.csv(Cheng.gene.all.anno.2.data.table.3,file = "/media/H_driver/2015/Nimer_Cheng/GeneWise_jscs3_all_with_anno_2_24_2016.csv")
#
# #GO term analysis
# #Adjust by the number of SJ
# Re.Go.adjusted.by.number.junction.22<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(data.table.gene.based.all.jscs4,ad="J","J",0.05,"Splice_junction_based")
# data.pwf2.SJs<-plotPWF2(Re.Go.adjusted.by.number.junction.22[[2]],binsize=30,xlab = "Number of SJs(<binsize> gene bins)")
# #Output annotated GO term
# capture.output(for(go in Re.Go.adjusted.by.number.junction.22[[3]]$category)
# { print(GOTERM[[go]])
#   cat("--------------------------------------\n")
# }, file="/media/H_driver/2015/Nimer_Cheng/SigGo_adjusted_by_SJ_2_24_2016.txt")
#
# Re.Go.adjusted.by.number.junction.sampling.based<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(data.table.gene.based.all.jscs4,ad="J","J",0.05,"Splice_junction_based")
# data.pwf2.SJs.sampling<-plotPWF2(Re.Go.adjusted.by.number.junction.sampling.based[[2]],binsize=30,xlab = "Number of SJs(<binsize> gene bins)")
# #Output annotated GO term
# capture.output(for(go in Re.Go.adjusted.by.number.junction.sampling.based[[3]]$category)
# { print(GOTERM[[go]])
#   cat("--------------------------------------\n")
# }, file="/media/H_driver/2015/Nimer_Cheng/SigGo_adjusted_by_SJ_sampling_based_2_25_2016.txt")
#
# #use the same data set from SJ analysis, but adjust analysis by gene length
# Re.Go.adjusted.by.number.junction.using.same.data.set<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(data.table.gene.based.all.jscs4,ad="GL","J",0.05,"Splice_junction_based")
# data.pwf2.adjust.by.gene.length.using.same.data.set<-plotPWF2(Re.Go.adjusted.by.number.junction.using.same.data.set[[2]],binsize=30,xlab = "gene length")
# #Output annotated GO term
# capture.output(for(go in Re.Go.adjusted.by.number.junction.using.same.data.set[[3]]$category)
# { print(GOTERM[[go]])
#   cat("--------------------------------------\n")
# }, file="/media/H_driver/2015/Nimer_Cheng/SigGo_adjusted_by_GL_using_same_data_set_as_SJ_analysis.txt")
#
# #Adjust by gene length
# Re.Go.adjusted.by.GL<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(data.table.gene.based.all.jscs4,ad="GL",NULL,0.05,"GL_based")
# data.pwf2.GLs<-plotPWF2(Re.Go.adjusted.by.GL[[2]],binsize=30,xlab = "gene length")
# #Output annotated GO term
# capture.output(for(go in Re.Go.adjusted.by.GL[[3]]$category)
# { print(GOTERM[[go]])
#   cat("--------------------------------------\n")
# }, file="/media/H_driver/2015/Nimer_Cheng/SigGo_adjusted_by_GL_2_25_2016.txt")
#
