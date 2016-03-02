#Install required packages
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("DESeq2")
biocLite("DEXSeq")
biocLite("JunctionSeq")
biocLite("KEGG.db")
library(KEGG.db)

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

#Run for the data from cheng
jscs.3<-GetResultsFromJunctionSeq("/media/H_driver/2015/Nimer_Cheng/","decoder.bySample.txt","_junction_seq_new_gtf_7/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt",
                          "/media/DATA/mus_musculus/Mus_musculus.GRCm38.83.JunctionSeq.flat.gff")
data.jscs.3<-fData(jscs.3)
gene.id.all.3<-unique(as.character(data.jscs.3[,1]))
table.gene.based.all.jscs3<-makeGeneWiseTable(jscs.3,gene.id.all.3,FDR.threshold = 0.05,verbose = TRUE, debug.mode = FALSE)
data.table.gene.based.all.jscs3<-pData(table.gene.based.all.jscs3)

#extract genes that are testable
data.table.gene.based.all.jscs4<-data.table.gene.based.all.jscs3[which(as.character(data.table.gene.based.all.jscs3[,8])!="character(0)"),]

#Annotate all genes
Cheng.gene.all.anno.new<-do.call(rbind,lapply(data.table.gene.based.all.jscs4[,1],GenerateGeneAnno))
Cheng.gene.all.anno.2.data.table.2<-merge(Cheng.gene.all.anno.3,data.table.gene.based.all.jscs4,by="geneID",sort=FALSE)

#Output results
Cheng.gene.all.anno.2.data.table.3<- data.frame(lapply(Cheng.gene.all.anno.2.data.table.2, as.character), stringsAsFactors=FALSE)
write.csv(Cheng.gene.all.anno.2.data.table.3,file = "/media/H_driver/2015/Nimer_Cheng/GeneWise_jscs3_all_with_anno_2_24_2016.csv")

#GO term analysis
#Adjust by the number of SJ
Re.Go.adjusted.by.number.junction.sampling.based<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(data.table.gene.based.all.jscs4,ad="J","J",0.05,"Splice_junction_based_2_25_2016")
data.pwf2.SJs.sampling<-plotPWF2(Re.Go.adjusted.by.number.junction.sampling.based[[2]],binsize=30,xlab = "Number of SJs(<binsize> gene bins)")

#Output annotated GO term
capture.output(for(go in Re.Go.adjusted.by.number.junction.sampling.based[[3]]$category)
{ print(GOTERM[[go]])
  cat("--------------------------------------\n")
}, file="/media/H_driver/2015/Nimer_Cheng/SigGo_adjusted_by_SJ_sampling_based_2_25_2016.txt")

#Adjust by the number of exon
Re.Go.adjusted.by.number.exon.sampling.based<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(data.table.gene.based.all.jscs4,ad="E","E",0.05,"Exon_based_2_26_2016")
data.pwf2.exons.sampling<-plotPWF2(Re.Go.adjusted.by.number.exon.sampling.based[[2]],binsize=30,xlab = "Number of exons(<binsize> gene bins)")

#Output annotated GO term
capture.output(for(go in Re.Go.adjusted.by.number.junction.sampling.based[[3]]$category)
{ print(GOTERM[[go]])
  cat("--------------------------------------\n")
}, file="/media/H_driver/2015/Nimer_Cheng/SigGo_adjusted_by_SJ_sampling_based_2_25_2016.txt")

#Adjust by the number of splice junction and exons
Re.Go.adjusted.by.number.exon.junction<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(data.table.gene.based.all.jscs4,ad="EJ",NULL,0.05,"Splice_exon_based_2_25_2016")
data.pwf2.exons.junction<-plotPWF2(Re.Go.adjusted.by.number.exon.junction[[2]],binsize=30,xlab = "Number of exons and splice junctions(<binsize> gene bins)")

#KEGG enrichment analysis
#adjust by the number of splice junction
Re.Go.adjusted.by.number.junction.sampling.based.kegg<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene_KEGG(data.table.gene.based.all.jscs4,ad="J","J",0.05,"Splice_junction_based_2_26_2016_2")
data.pwf2.SJs.sampling<-plotPWF2(Re.Go.adjusted.by.number.junction.sampling.based.kegg[[2]],binsize=30,xlab = "Number of SJs(<binsize> gene bins)")

#adjust by the number of exons
Re.Go.adjusted.by.number.exon.sampling.based.kegg<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene_KEGG(data.table.gene.based.all.jscs4,ad="E","E",0.05,"Exon_based_2_26_2016")
data.pwf2.exons.sampling<-plotPWF2(Re.Go.adjusted.by.number.exon.sampling.based[[2]],binsize=30,xlab = "Number of exons(<binsize> gene bins)")

#adjust by the number of exons and splice junctions
Re.Go.adjusted.by.number.exon.junction<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene_KEGG(data.table.gene.based.all.jscs4,ad="EJ",NULL,0.05,"Splice_exon_based_2_26_2016")
data.pwf2.exons.junction<-plotPWF2(Re.Go.adjusted.by.number.exon.junction[[2]],binsize=30,xlab = "Number of exons and splice junctions(<binsize> gene bins)")

#Anotate KEGG pathway
keggid2keggname <- as.list(KEGGPATHID2NAME)
AddKEGGAnnotation<-function(Data2BeAnno,keggName,file_prefix){

keggid2<-Data2BeAnno[[1]][,1]
Data2BeAnno.annotation<-as.data.frame(unlist(keggid2keggname[keggid2]))
Data2BeAnno.annotation.2<-cbind(rownames(Data2BeAnno.annotation),as.character(Data2BeAnno.annotation[,1]))
colnames(Data2BeAnno.annotation.2)=c("category","annotation")
head(Data2BeAnno.annotation.2)
head(Data2BeAnno[[1]])
Data2BeAnno.anno<-merge(Data2BeAnno[[1]],Data2BeAnno.annotation.2,by=1)

keggid3<-Data2BeAnno[[3]][,1]
Data2BeAnno.annotation.enriched<-as.data.frame(unlist(keggid2keggname[keggid3]))
Data2BeAnno.annotation.enriched.2<-cbind(rownames(Data2BeAnno.annotation.enriched),as.character(Data2BeAnno.annotation.enriched[,1]))
colnames(Data2BeAnno.annotation.enriched.2)=c("category","annotation")
head(Data2BeAnno.annotation.2)
head(Data2BeAnno[[1]])
Data2BeAnno.anno.enriched<-merge(Data2BeAnno[[3]],Data2BeAnno.annotation.enriched.2,by=1)

write.csv(Data2BeAnno.anno,file=paste(file_prefix,"_all_KEGG",".csv",sep=""),row.names = FALSE)
write.csv(Data2BeAnno.anno.enriched,file=paste(file_prefix,"_enriched_KEGG",".csv",sep=""),row.names = FALSE)

}

#Add annotations to the KEGG enrichment analysis adjusted by the number of splice junctions
AddKEGGAnnotation(Re.Go.adjusted.by.number.junction.sampling.based.kegg,keggid2keggname,"Splice_junction_based")

#Add annotations to the KEGG enrichment analysis adjusted by the number of exons
AddKEGGAnnotation(Re.Go.adjusted.by.number.exon.sampling.based.kegg,keggid2keggname,"Exon_junction_based")

#Add annotations to the KEGG enrichment analysis adjusted by the number of splice junctions and exons
AddKEGGAnnotation(Re.Go.adjusted.by.number.exon.junction,keggid2keggname,"Splice_junction_exon_based")

#Re.Go.adjusted.by.number.exon.junction[[2]]
#goseq.m(Re.Go.adjusted.by.number.exon.junction[[2]],"mm10","ensGene",method = "Sampling",repcnt=10000,test.cats = c("KEGG"))

#test.KEGG<-as.data.frame(unlist(getgo.m(rownames(Re.Go.adjusted.by.number.exon.junction[[2]]),"mm10","ensGene",fetch.cats=c("KEGG"))))
#test.KEGG<-getgo.m(rownames(Re.Go.adjusted.by.number.exon.junction[[2]]),"mm10","ensGene",fetch.cats=c("KEGG"))
