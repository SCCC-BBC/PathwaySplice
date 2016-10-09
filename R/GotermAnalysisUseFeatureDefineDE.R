#' GotermAnalysisUseFeatureDefineDE
#' Descrption: This function computes statistical significance of GO terms, given a list of significant genes. These significant genes can be directly from differntial splicing analysis software such as rMat, DEXSeq or JunctionSeq. 

#' Based on genewise results to perform Go term analysis using DGEs defined by different criterion
#' @param Data4Goterm
#' @param ad
#' @param sub_feature
#' @param gene_model
#'
#'
#' @return
#' @export
#'
#' @examples
#'
#'example.GO<-GotermAnalysisUseFeatureDefineDE(mds,ad="exon_SJ",sub_feature="E",DE_define="Feature",gene_model="hg19",Output_file_dir=paste0(getwd(),"/"))
#'
GotermAnalysisUseFeatureDefineDE<-function(re.gene.based,ad="GL",sub_feature=NULL,DE_define=c("Feature","GeneWise","rMAT","FeatureGeneWise","FeaturerMAT","GeneWiserMAT","FeatureGeneWiseRMAT"),gene_model,Output_file_dir){

  #Data4Goterm<-pData(re.gene.based)

  Data4Goterm<-re.gene.based

  if(is.null(sub_feature)){Data4Goterm.sub_feature<-Data4Goterm}
  else{Data4Goterm.sub_feature<-Data4Goterm[grep(sub_feature,Data4Goterm[,8]),]}

  #print(dim(Data4Goterm.sub_feature))

  if(sub_feature=="J"){
  Data4Goterm.sub_feature.geneID.NumOfJunctions<-Data4Goterm.sub_feature[,c(1,11)]
  }else
  {
    Data4Goterm.sub_feature.geneID.NumOfJunctions<-Data4Goterm.sub_feature[,c(1,10)]
  }
  #print(dim(Data4Goterm.sub_feature.geneID.NumOfJunctions))

  #Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature[,7]<threshold),]

  if(DE_define=="Feature"){
  Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature$DE_or_not_feature==1),]
  }else if(DE_define=="GeneWise"){
  Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature$DE_or_not_geneWise==1),]
  }else if(DE_define=="rMAT"){
    Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature$DE_or_not_rMAT_based==1),]
  }else if(DE_define=="FeatureGeneWise"){
  Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature$DE_or_not_feature==1&Data4Goterm.sub_feature$DE_or_not_geneWise==1),]
  }else if(DE_define=="FeaturerMAT"){
    Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature$DE_or_not_feature==1&
                                                                 Data4Goterm.sub_feature$DE_or_not_rMAT_based==1),]
  }else if(DE_define=="GeneWiserMAT"){
    Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature$DE_or_not_geneWise==1&
                                                               Data4Goterm.sub_feature$DE_or_not_rMAT_based==1),]
  }else if(DE_define=="FeatureGeneWiseRMAT"){
  Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature$DE_or_not_feature==1&
                                                               Data4Goterm.sub_feature$DE_or_not_geneWise==1&
                                                               Data4Goterm.sub_feature$DE_or_not_rMAT_based==1),]}

  #GO term analysis using GOSeq

  All.gene.id.based.on.sub_feature<-unique(Data4Goterm.sub_feature[,1])

  cat("How many DE genes?","\n")
  cat(dim(Data4Goterm.sub_feature.Sig)[1],"\n")
  cat(length(unique(as.character(Data4Goterm.sub_feature.Sig$geneID))),"\n")

  All.gene.id.index<-rep(0,length(All.gene.id.based.on.sub_feature))
  names(All.gene.id.index)=All.gene.id.based.on.sub_feature

  All.genes.based.on.Sig.sub_feature<-unique(Data4Goterm.sub_feature.Sig[,1])
  gene.DE_interest<-as.integer(which( All.gene.id.based.on.sub_feature %in% All.genes.based.on.Sig.sub_feature ))

  All.gene.id.index[gene.DE_interest]<-1
  #print(length(All.gene.id.index))

  gene.with.matched.junction<-which(Data4Goterm.sub_feature.geneID.NumOfJunctions[,1] %in% c(names(All.gene.id.index)))
  num.junction.4.matched.gene<-as.numeric(Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,2])

  #names.4.matched.gene<-Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,1]

  #All.gene.id.index.2<-All.gene.id.index[which(names(All.gene.id.index) %in% c(names.4.matched.gene))]

  #print(length(All.gene.id.index.2))

  All.gene.id.index.2<-All.gene.id.index

  #print(All.gene.id.index.2)

  if(ad=="GL"){
    pwf.DE_interest=nullp(All.gene.id.index.2,gene_model,"ensGene",plot.fit = FALSE)
  }
  else
  {
    pwf.DE_interest=nullp(All.gene.id.index.2,gene_model,"ensGene",bias.data = num.junction.4.matched.gene,plot.fit = FALSE)
  }

  GO.wall.DE_interest=goseq2(pwf.DE_interest,gene_model,"ensGene",gene_model,test.cats=c("GO:BP"),use_genes_without_cat=TRUE)

  #GO.wall.DE_interest=goseq(pwf.DE_interest,"mm10","ensGene",test.cats=c("GO:BP"),use_genes_without_cat=TRUE)

  #GO.wall.DE_interest=goseq2(pwf.DE_interest,"mm10","ensGene",gene.model=gene_model)
  #enriched.GO.DE_interest=GO.wall.DE_interest[p.adjust(GO.wall.DE_interest$over_represented_pvalue,method="BH")<threshold,]
  re<-list()

  re[[1]]<-GO.wall.DE_interest[[1]]
  re[[2]]<-pwf.DE_interest
  re[[3]]<-GO.wall.DE_interest[[2]]


  DE.gene.symbol<-gene.model[which(as.character(gene.model[,3]) %in% unique(as.character(Data4Goterm.sub_feature.Sig$geneID))),1]

  Re4_temp<-list(A=DE.gene.symbol,B=re[[3]])

  venn.plot <- venn.diagram(
    x = Re4_temp[c(1,2)],
    filename = paste0(Output_file_dir,DE_define,"_overlap_with_DE_from_GO.tiff"),
    col = "black",
    lty = "dotted",
    lwd = 2,
    fill = c("red","blue"),
    alpha = 0.50,
    label.col = c(rep("white",3)),
    cex = 1,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("red","blue"),
    cat.cex = 0.8,
    cat.fontfamily = "serif"
  )

  return(re)
}
