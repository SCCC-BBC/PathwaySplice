#' GetPwfUseReformatedData
#'
#' Use reformated data to calculate probability weight function
#'
#' @param re.gene.based
#' @param ad
#' @param sub_feature
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
#' Re.pwf.exon.sj<-GetPwfUseReformatedData(mds,ad="exon_SJ",sub_feature="E",0.05)
#' 
GetPwfUseReformatedData<-function(re.gene.based,ad="GL",sub_feature=NULL,threshold){

  Data4Goterm<-re.gene.based

  if(is.null(sub_feature)){Data4Goterm.sub_feature<-Data4Goterm}
  else{Data4Goterm.sub_feature<-Data4Goterm[grep(sub_feature,Data4Goterm[,8]),]}

  Data4Goterm.sub_feature.geneID.NumOfJunctions<-Data4Goterm.sub_feature[,c(1,10)]

  Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature[,7]<threshold),]

  All.gene.id.based.on.sub_feature<-unique(Data4Goterm.sub_feature[,1])
  All.gene.id.index<-rep(0,length(All.gene.id.based.on.sub_feature))
  names(All.gene.id.index)=All.gene.id.based.on.sub_feature

  All.genes.based.on.Sig.sub_feature<-unique(Data4Goterm.sub_feature.Sig[,1])
  gene.DE_interest<-as.integer(which( All.gene.id.based.on.sub_feature %in% All.genes.based.on.Sig.sub_feature ))

  All.gene.id.index[gene.DE_interest]<-1
  #print(length(All.gene.id.index))

  gene.with.matched.junction<-which(Data4Goterm.sub_feature.geneID.NumOfJunctions[,1] %in% c(names(All.gene.id.index)))
  num.junction.4.matched.gene<-as.numeric(Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,2])

  All.gene.id.index.2<-All.gene.id.index

  #print(All.gene.id.index.2)

  if(ad=="GL"){
    pwf.DE_interest=nullp(All.gene.id.index.2,"mm10","ensGene",plot.fit = FALSE)
  }
  else
  {
    pwf.DE_interest=nullp(All.gene.id.index.2,"mm10","ensGene",bias.data = num.junction.4.matched.gene,plot.fit = FALSE)
  }

  re<-pwf.DE_interest

  return(re)
}
