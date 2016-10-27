#' Run_pathwaysplice
#'
#' @param Data4Goterm
#' @param ad
#' @param sub_feature
#' @param threshold
#'
#' @import goseq
#'
#' @importFrom geneLenDataBase unfactor
#' @importFrom BiasedUrn dWNCHypergeo pWNCHypergeo
#' @importFrom gdata trim
#' @importFrom GO.db GO.db
#'
#' @export
#'
#' @examples
#' data(mds)
#' data(hg19)
#'
#' Example.Go.adjusted.by.exon<-Run_pathwaysplice(mds,ad="exon_SJ",sub_feature="E",0.05,genomeID="hg19",geneID="ensGene",gene_model=hg19.gene.model,method="Wallenius")
#'
#' Example.Go.unadjusted<-Run_pathwaysplice(mds,ad="exon_SJ",sub_feature="E",0.05,genomeID="hg19",geneID="ensGene",gene_model=hg19.gene.model,method="Hypergeometric")
#'
#' write.table(Example.Go.unadjusted[[1]]$DE_GO,file="DE.txt",col.names = F,row.names = F,quote=F)
#'
Run_pathwaysplice<-function(re.gene.based,ad="GL",sub_feature=NULL,threshold,genomeID,geneID,gene_model,method){

  #Data4Goterm<-pData(re.gene.based)

  Data4Goterm<-re.gene.based

  if(is.null(sub_feature)){Data4Goterm.sub_feature<-Data4Goterm}
  else{Data4Goterm.sub_feature<-Data4Goterm[grep(sub_feature,Data4Goterm[,8]),]}


  #print(dim(Data4Goterm.sub_feature))

  if(sub_feature=="J"){
  Data4Goterm.sub_feature.geneID.NumOfJunctions<-Data4Goterm.sub_feature[,c(1,11)]
  }else{
    Data4Goterm.sub_feature.geneID.NumOfJunctions<-Data4Goterm.sub_feature[,c(1,10)]
  }
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

  if(method=="Hypergeometric"){
  GO.wall.DE_interest=pathwaysplice(pwf.DE_interest,genomeID,geneID,gene.model=gene_model,method="Hypergeometric",use_genes_without_cat=TRUE)
  }else
  {
  GO.wall.DE_interest=pathwaysplice(pwf.DE_interest,genomeID,geneID,gene.model=gene_model,use_genes_without_cat=TRUE)
  }
  #GO.wall.DE_interest=goseq2(pwf.DE_interest,"mm10","ensGene",gene.model=gene_model)
  #enriched.GO.DE_interest=GO.wall.DE_interest[p.adjust(GO.wall.DE_interest$over_represented_pvalue,method="BH")<threshold,]
  re<-list(GO.wall.DE_interest=GO.wall.DE_interest,pwf.DE_interest)

  #re[[1]]<-GO.wall.DE_interest
  #re[[2]]<-pwf.DE_interest

  return(re)
}
