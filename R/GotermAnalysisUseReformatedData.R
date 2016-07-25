#' GotermAnalysisUseReformatedData
#'
#'
#' @param Data4Goterm
#' @param ad
#' @param sub_feature
#' @param threshold
#'
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' Re.Go.adjusted.by.number.junction.2<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(
#' re.PJ.gene.based,ad="J","J",0.05,"Splice_junction_based")
#'
#' data.pwf2.SJs<-plotPWF2(Re.Go.adjusted.by.number.junction.2[[2]],binsize=30,xlab = "Number of SJs(<binsize> gene bins)")
#'
#' Re.Go.adjusted.by.exon.SJ<-GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene(
#' re.PJ.gene.based,ad="exon_SJ",sub_feature=NULL,0.05,file_prefix="Exon_Splice_junction_based.xls",gene_model=gene.model)
#'
GotermAnalysisUseReformatedData<-function(re.gene.based,ad="GL",sub_feature=NULL,threshold,genomeID,geneID,gene_model){

  #Data4Goterm<-pData(re.gene.based)

  Data4Goterm<-re.gene.based

  if(is.null(sub_feature)){Data4Goterm.sub_feature<-Data4Goterm}
  else{Data4Goterm.sub_feature<-Data4Goterm[grep(sub_feature,Data4Goterm[,8]),]}


  #print(dim(Data4Goterm.sub_feature))

  Data4Goterm.sub_feature.geneID.NumOfJunctions<-Data4Goterm.sub_feature[,c(1,11)]

  #print(dim(Data4Goterm.sub_feature.geneID.NumOfJunctions))

  Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature[,7]<threshold),]

  #GO term analysis using GOSeq
  All.gene.id.based.on.sub_feature<-unique(Data4Goterm.sub_feature[,1])
  #length(All.gene.id.based.on.sub_feature)
  All.gene.id.index<-rep(0,length(All.gene.id.based.on.sub_feature))
  names(All.gene.id.index)=All.gene.id.based.on.sub_feature

  All.genes.based.on.Sig.sub_feature<-unique(Data4Goterm.sub_feature.Sig[,1])
  gene.DE_interest<-as.integer(which( All.gene.id.based.on.sub_feature %in% All.genes.based.on.Sig.sub_feature ))

  All.gene.id.index[gene.DE_interest]<-1
  print(length(All.gene.id.index))

  gene.with.matched.junction<-which(Data4Goterm.sub_feature.geneID.NumOfJunctions[,1] %in% c(names(All.gene.id.index)))
  num.junction.4.matched.gene<-as.numeric(Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,2])

  #names.4.matched.gene<-Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,1]

  #All.gene.id.index.2<-All.gene.id.index[which(names(All.gene.id.index) %in% c(names.4.matched.gene))]

  #print(length(All.gene.id.index.2))

  All.gene.id.index.2<-All.gene.id.index

  print(All.gene.id.index.2)

  if(ad=="GL"){
    pwf.DE_interest=nullp(All.gene.id.index.2,genomeID,geneID,plot.fit = FALSE)
  }
  else
  {
    pwf.DE_interest=nullp(All.gene.id.index.2,genomeID,geneID,bias.data = num.junction.4.matched.gene,plot.fit = FALSE)
  }

  GO.wall.DE_interest=goseq2(pwf.DE_interest,genomeID,geneID,gene.model=gene_model,use_genes_without_cat=TRUE)

  #GO.wall.DE_interest=goseq2(pwf.DE_interest,"mm10","ensGene",gene.model=gene_model)
  enriched.GO.DE_interest=GO.wall.DE_interest[p.adjust(GO.wall.DE_interest$over_represented_pvalue,method="BH")<threshold,]
  re<-list()

  re[[1]]<-GO.wall.DE_interest
  re[[2]]<-pwf.DE_interest

  return(re)
}
