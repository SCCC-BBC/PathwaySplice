#' Select_DE_gene_based_on_Feature
#'
#' @param Re.PJ
#' @param re.PJ.gene.based
#' @param re.rMAT
#' @param splicing_type: there are 5 splicing types from rMAT, need to choose one
#'
#' @return
#' @export
#'
#' @examples
#'
#' Re.PJ.selected.feature.FC.p<-Select_DE_gene_basd_on_Feature(Re.PJ,re.PJ.gene.based,re.rMAT,0.58,0.05,outputfile_DGE_FC_P_geneWise)
#' dim(Re.PJ.selected.feature.FC.p)
#' head(Re.PJ.selected.feature.FC.p)
#' length(which(Re.PJ.selected.feature.FC.p$DE_or_not==1))
#'

Select_DE_gene_basd_on_Feature<-function(Re.PJ,re.PJ.gene.based,re.rMAT,splicing_type,cutoff_FC,cutoff_P_value,venn_output){

  feature.based.RE.PJ<-fData(Re.PJ)

  gene.based.RE.PJ<-ReformatDataUseTable(pData(re.PJ.gene.based))

  cat(dim(feature.based.RE.PJ),"\t",dim(gene.based.RE.PJ),"\n")

  re2<-feature.based.RE.PJ[which(feature.based.RE.PJ[,20]>cutoff_FC&feature.based.RE.PJ[,11]<cutoff_P_value),]

  DE.gene.based.on.FC.p.of.feature<-unique(unlist(strsplit(as.character(re2$geneID),"\\+")))

  re2.geneWise.p.based<-feature.based.RE.PJ[which(feature.based.RE.PJ$geneWisePadj<cutoff_P_value),]

  DE.gene.based.on.geneWise.p.only<-unique(unlist(strsplit(as.character(re2.geneWise.p.based$geneID),"\\+")))

  #DE.gene.based.on.FC.p.of.feature.and.geneWise.p<-intersect(DE.gene.based.on.FC.p.of.feature,DE.gene.based.on.geneWise.p.only)



  DE.gene.feature<-DE.gene.based.on.FC.p.of.feature
  DE.gene.geneWisePadj<-DE.gene.based.on.geneWise.p.only

  if(splicing_type=="All_5_Types"){
    cat(length(DE.gene.based.on.FC.p.of.feature),"\t",length(DE.gene.based.on.geneWise.p.only),"\t",
        length(re.rMAT$ReadsOnTargetAndJunctionCounts),"\n")
  DE.gene.rMAT<-re.rMAT$ReadsOnTargetAndJunctionCounts
  }else if(splicing_type=="SE"){

    cat(length(DE.gene.based.on.FC.p.of.feature),"\t",length(DE.gene.based.on.geneWise.p.only),"\t",
        length(re.rMAT$SEReadsOnTargetAndJunctionCounts),"\n")
      DE.gene.rMAT<-re.rMAT$SEReadsOnTargetAndJunctionCounts
    }

  DE_or_not_feature<-rep(0,dim(gene.based.RE.PJ)[1])
  DE_or_not_geneWise<-rep(0,dim(gene.based.RE.PJ)[1])
  DE_or_not_rMAT_based<-rep(0,dim(gene.based.RE.PJ)[1])

  re3<-cbind(gene.based.RE.PJ,DE_or_not_feature,DE_or_not_geneWise,DE_or_not_rMAT_based)


  re3[which(re3$geneID %in% DE.gene.feature),]$DE_or_not_feature<-1
  re3[which(re3$geneID %in% DE.gene.geneWisePadj),]$DE_or_not_geneWise<-1
  re3[which(re3$geneID %in% DE.gene.rMAT),]$DE_or_not_rMAT_based<-1

  #re4<-unique(re3[which(re3$DE_or_not==1),]$geneID)

  Re4<-list(GeneAll=re3,DGE.geneWise=DE.gene.geneWisePadj,
            DGE.featureWise=DE.gene.feature,
            DGE.rMAT=DE.gene.rMAT)

  # venn.plot <- venn.diagram(
  #   x = Re4[c(1,2)],
  #   filename = venn_output,
  #   height = 3000,
  #   width = 3500,
  #   resolution = 1000,
  #   col = "black",
  #   lty = "dotted",
  #   lwd = 1,
  #   fill = c("red","blue"),
  #   alpha = 0.50,
  #   label.col = c(rep("white",3)),
  #   cex = 0.5,
  #   fontfamily = "serif",
  #   fontface = "bold",
  #   cat.col = c("red","blue"),
  #   cat.cex = 0.5,
  #   cat.pos = 0.5,
  #   cat.dist = 0.05,
  #   cat.fontfamily = "serif"
  # )

  venn.plot <- venn.diagram(
    x = Re4[c(2,3,4)],
    filename = paste0(venn_output,"DE_gene_overlap_rMAT_SE.tiff"),
    col = "black",
    lty = "dotted",
    lwd = 2,
    fill = c("red", "orange", "blue"),
    alpha = 0.50,
    label.col = c(rep("white",7)),
    cex = 1,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("red", "orange", "blue"),
    cat.cex = 0.8,
    cat.fontfamily = "serif"
  )

  write.table(as.data.frame(intersect(Re4$DGE.featureWise,Re4$DGE.rMAT)),file=paste0(venn_output,"DGE_overlap_rMAT_SE_with_feature.txt"),row.names = FALSE,quote=FALSE,sep="\t")

  return(Re4)

}
