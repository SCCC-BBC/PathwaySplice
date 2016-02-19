GotermAnalysis2GetAllGOTerms_AdjustedByNumOfJunctionWithinOneGene<-function(Data4Goterm,ad="GL",sub_feature=NULL,threshold,file_prefix){
  
  if(is.null(sub_feature)){Data4Goterm.sub_feature<-Data4Goterm}
  else{Data4Goterm.sub_feature<-Data4Goterm[grep(sub_feature,Data4Goterm[,8]),]}

    
  print(dim(Data4Goterm.sub_feature))
  
  Data4Goterm.sub_feature.geneID.NumOfJunctions<-Data4Goterm.sub_feature[,c(1,11)]
  
  print(dim(Data4Goterm.sub_feature.geneID.NumOfJunctions))
    
  Data4Goterm.sub_feature.Sig<-Data4Goterm.sub_feature[which(Data4Goterm.sub_feature[,7]<threshold),]
  
  #dim(Data4Goterm.sub_feature.Sig)
  
  #unique(Data4Goterm.sub_feature.Sig[,1])
  
  #GO term analysis using GOSeq
  All.gene.id.based.on.sub_feature<-unique(unlist(strsplit(Data4Goterm.sub_feature[,1],"\\+")))
  #length(All.gene.id.based.on.sub_feature)
  All.gene.id.index<-rep(0,length(All.gene.id.based.on.sub_feature))
  names(All.gene.id.index)=All.gene.id.based.on.sub_feature
  
  All.genes.based.on.Sig.sub_feature<-unique(unlist(strsplit(Data4Goterm.sub_feature.Sig[,1],"\\+")))
  #length(All.genes.based.on.Sig.sub_feature)
  
  #gene.de<-intersect(All.gene.id.based.on.sub_feature,All.genes.based.on.Sig.sub_feature)
  #length(gene.de)
  
  gene.DE_interest<-as.integer(which( All.gene.id.based.on.sub_feature %in% All.genes.based.on.Sig.sub_feature ))
  
  #gene.DE_interest
  #All.gene.id.index
  All.gene.id.index[gene.DE_interest]<-1
  #All.gene.id.index
  
  #names(All.gene.id.index)
  #names(All.gene.id.index)<-All.gene.id.based.on.sub_feature
  
  #genes2go.DE_interest=getgo(names(All.gene.id.index),'mm10','ensGene')
  #genes2go.DE_interest
  #go2genes.DE_interest=goseq:::reversemapping(genes2go.DE_interest) 
  #go2genes.DE_interest
  #print(names(All.gene.id.index))
  
  gene.with.matched.junction<-which(Data4Goterm.sub_feature.geneID.NumOfJunctions[,1] %in% c(names(All.gene.id.index)))
  
  num.junction.4.matched.gene<-Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,2]
  
  names.4.matched.gene<-Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction,1]
  
  All.gene.id.index.2<-All.gene.id.index[which(names(All.gene.id.index) %in% c(names.4.matched.gene))]
  
#  print(dim(All.gene.id.index))
  print(length(All.gene.id.index.2))   
  
  
if(ad=="GL"){
pwf.DE_interest=nullp(All.gene.id.index.2,"mm10","ensGene",plot.fit = FALSE)
}
else
{
pwf.DE_interest=nullp(All.gene.id.index.2,"mm10","ensGene",bias.data = num.junction.4.matched.gene,plot.fit = FALSE)  
}
  
  #head(pwf.DE_interest)
  #dim(pwf.DE_interest)
  #length(which(pwf.DE_interest[,1]==1))
  #length(which(pwf.DE_interest[,1]==0))
  
  GO.wall.DE_interest=goseq(pwf.DE_interest,"mm10","ensGene")
  #dim(GO.wall.DE_interest)
  #head(GO.wall.DE_interest)
  enriched.GO.DE_interest=GO.wall.DE_interest[p.adjust(GO.wall.DE_interest$over_represented_pvalue,method="BH")<threshold,]
  
  #  enriched.GO.DE_interest=GO.wall.DE_interest$category[p.adjust(GO.wall.DE_interest$over_represented_pvalue,method="BH")<threshold]
  
  #  capture.output(for(go in enriched.GO.DE_interest) { print(GOTERM[[go]])
  #    cat("--------------------------------------\n")
  #  }, file=paste(file_prefix,"_",threshold,"_SigGo",".txt",sep=""))
  
  # enriched.GO.DE_interest.005=GO.wall.DE_interest$category[p.adjust(GO.wall.DE_interest$over_represented_pvalue,method="BH")<0.05]
  # 
  # capture.output(for(go in enriched.GO.DE_interest.005) { print(GOTERM[[go]])
  #   cat("--------------------------------------\n")
  # }, file="SigGo_based_on_splice_juncton_0.05.txt")
  
 #write.csv(GO.wall.DE_interest,file=paste(file_prefix,"_",threshold,"_Go_adjusted_by_num_of_junctions",".csv",sep=""),row.names = FALSE)
  
  re<-list()
  
  re[[1]]<-GO.wall.DE_interest
  re[[2]]<-pwf.DE_interest
  return(re)

  }
