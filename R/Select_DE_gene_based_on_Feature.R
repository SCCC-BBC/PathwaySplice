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
#' Re.PJ.selected.feature.FC.p.check<-Select_DE_gene_basd_on_Feature(Re.PJ,re.PJ.gene.based,re.rMAT,"SE",2,0.05,outputfile_DGE_FC_P_geneWise)

Select_DE_gene_basd_on_Feature<-function(Re.PJ,re.PJ.gene.based,re.rMAT,splicing_type,cutoff_FC,cutoff_P_value,gene_model,venn_output){

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

  FindElementsBetweenSets<-function(xx,gene_model,venn_output){

    a1<-xx[[1]]
    a2<-xx[[2]]
    a3<-xx[[3]]

    a12<-intersect(a1,a2)
    a13<-intersect(a1,a3)
    a23<-intersect(a2,a3)

    a123<-intersect(a12,a3)
    a12not3<-setdiff(a12,a123)
    a13not2<-setdiff(a13,a123)
    a23not1<-setdiff(a23,a123)

    a1not23<-setdiff(a1,union(a12,a13))
    a2not13<-setdiff(a2,union(a12,a23))
    a3not12<-setdiff(a3,union(a13,a23))

    DGE.geneWise<-c(a1not23,a12not3,a13not2,a123)
    DGE.featureWise<-c(a12not3,a2not13,a123,a23not1)
    DGE.rMAT<- c(a13not2,a123,a23not1,a3not12)

    re<-qpcR:::cbind.na(a123,a12not3,a13not2,a23not1,a1not23,a2not13,a3not12,DGE.geneWise,DGE.featureWise,DGE.rMAT)
    re[is.na(re)] <- ""
    re<-data.frame(apply(re,2,sort,decreasing=T))

    colnames(re)<-c(paste0(length(a123),"_genes"),paste0(length(a12not3),"_genes"),paste0(length(a13not2),"_genes"),
                    paste0(length(a23not1),"_genes"),paste0(length(a1not23),"_genes"),paste0(length(a2not13),"_genes"),
                    paste0(length(a3not12),"_genes"),paste0(length(DGE.geneWise),"_DGE.geneWise"),
                    paste0(length(DGE.featureWise),"_DGE.featureWise"),paste0(length(DGE.rMAT),"_DGE.rMAT"))

     a123.gene.symbol<-sapply(a123,function(u,gene.model){
      y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
      if(length(as.character(y))==0){
        y<-"NotMatched"
      }
      y
    },gene.model)

     a12not3.gene.symbol<-sapply(a12not3,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)

     a13not2.gene.symbol<-sapply(a13not2,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)

     a23not1.gene.symbol<-sapply(a23not1,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)

     a1not23.gene.symbol<-sapply(a1not23,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)

     a2not13.gene.symbol<-sapply(a2not13,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)

     a3not12.gene.symbol<-sapply(a3not12,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)


     DGE.geneWise.gene.symbol<-sapply(DGE.geneWise,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)

     DGE.featureWise.gene.symbol<-sapply(DGE.featureWise,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)

     DGE.rMAT.gene.symbol<-sapply(DGE.rMAT,function(u,gene.model){
       y<-gene.model[which(as.character(gene.model[,3]) %in% u),1]
       if(length(as.character(y))==0){
         y<-"NotMatched"
       }
       y
     },gene.model)


     re2<-qpcR:::cbind.na(a123.gene.symbol,a12not3.gene.symbol,a13not2.gene.symbol,
                a23not1.gene.symbol,a1not23.gene.symbol,a2not13.gene.symbol,a3not12.gene.symbol,DGE.geneWise.gene.symbol,
                DGE.featureWise.gene.symbol,DGE.rMAT.gene.symbol)
     re2[is.na(re2)] <- ""
     re2<-data.frame(apply(re2,2,sort,decreasing=T))

     colnames(re2)<-c(paste0(length(a123.gene.symbol),"_genes"),paste0(length(a12not3.gene.symbol),"_genes"),paste0(length(a13not2.gene.symbol),"_genes"),
                     paste0(length(a23not1.gene.symbol),"_genes"),paste0(length(a1not23.gene.symbol),"_genes"),paste0(length(a2not13.gene.symbol),"_genes"),
                     paste0(length(a3not12.gene.symbol),"_genes"),
                     paste0(length(DGE.geneWise),"_DGE.geneWise"),
                     paste0(length(DGE.featureWise),"_DGE.featureWise"),paste0(length(DGE.rMAT),"_DGE.rMAT"))


    write.table(re,file=paste0(venn_output,"DGE_list_Ensemble.xls"),row.names = FALSE,quote=FALSE,sep="\t")

    write.table(re2,file=paste0(venn_output,"DGE_list_gene_symbol.xls"),row.names = FALSE,quote=FALSE,sep="\t")

   re3<-list(re=re,re2=re2)

   return(re3)

  }

  Ree4<-FindElementsBetweenSets(Re4[c(2,3,4)],gene_model,venn_output)


  return(Re4)

}
