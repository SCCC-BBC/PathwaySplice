library(PathwaySplice)
#Venn

PostProcessGO(25,Example.cp.adjusted.by.exon,Example.cp.unadjusted,
"/Volumes/Bioinformatics$/Aimin_project/ToGaoZhen/","/Volumes/Bioinformatics$/Aimin_project/ToGaoZhen/","In_ad_not_un.xls","In_un_not_ad.xls")


# n=25
# 
# Example.Go.adjusted.by.exon<-Example.cp.adjusted.by.exon
# Example.Go.unadjusted<-Example.cp.unadjusted
# 
# #Example.cp.adjusted.by.exon
# 
# adjusted<-Example.Go.adjusted.by.exon$GO.selected[1:n,1]
# unadjusted<-Example.Go.unadjusted$GO.selected[1:n,1]
# 
# re<-list(adjusted=adjusted,unadjusted=unadjusted)
# 
# venn.plot <- venn.diagram(
#   x = re[c(1,2)],
#   filename = paste0(getwd(),"/",names(re)[1],"_",names(re)[2],"_overlap_venn.tiff"),
#   #filename=NULL,
#   height = 3000,
#   width = 3500,
#   resolution = 1000,
#   col = "black",
#   lty = "dotted",
#   lwd = 1,
#   fill = c("red","blue"),
#   alpha = 0.50,
#   label.col = c(rep("black",3)),
#   cex = 0.5,
#   fontfamily = "serif",
#   fontface = "bold",
#   cat.col = c("red","blue"),
#   cat.cex = 0.5,
#   cat.pos = 0.5,
#   cat.dist = 0.05,
#   cat.fontfamily = "serif"
# )
# 
# 
# #boxplot
# 
# common<-intersect(unadjusted,adjusted)
# 
# In.unadjusted.not.in.adjusted<-setdiff(unadjusted,common)
# In.adjusted.not.in.unadjusted<-setdiff(adjusted,common)
# 
# length(In.unadjusted.not.in.adjusted)
# length(In.adjusted.not.in.unadjusted)
# length(common)
# 
# index1<-match(In.adjusted.not.in.unadjusted,Example.Go.adjusted.by.exon$GO.selected$category)
# In.ad.not.un<-Example.Go.adjusted.by.exon$GO.selected[index1,]$Ave_value_all_gene
# 
# index2<-match(In.unadjusted.not.in.adjusted,Example.Go.unadjusted$GO.selected$category)
# In.un.not.ad<-Example.Go.unadjusted$GO.selected[index2,]$Ave_value_all_gene
# 
# xx<-cbind(unlist(In.ad.not.un),unlist(In.un.not.ad))
# 
# colnames(xx)=c("In.ad.not.un","In.un.not.ad")
# 
# #boxplot(xx)
# 
# #cbind(Example.Go.adjusted.by.exon$GO.selected[index1,1],Example.Go.unadjusted$GO.selected[index2,1])
# 
# #In.ad.not.un<-xx[,1]
# #In.un.not.ad<-xx[,2]
# 
# cp.top.adjusted.25<-unlist(Example.Go.adjusted.by.exon$GO.selected[1:n,]$Ave_value_all_gene)
# cp.top.unadjusted.25<-unlist(Example.Go.unadjusted$GO.selected[1:n,]$Ave_value_all_gene)
# 
# cp.all.adjusted<-unlist(Example.Go.adjusted.by.exon$GO.selected$Ave_value_all_gene)
# cp.all.unadjusted<-unlist(Example.Go.unadjusted$GO.selected$Ave_value_all_gene)
# 
# yy<-rbind(cbind(xx[,1],rep("In.ad.not.un",length(xx[,1]))),
# cbind(xx[,2],rep("In.un.not.ad",length(xx[,2]))),
# cbind(cp.top.adjusted.25,rep("cp.top.adjusted.25",length(cp.top.adjusted.25))),
# cbind(cp.top.unadjusted.25,rep("cp.top.unadjusted.25",length(cp.top.unadjusted.25))),
# cbind(cp.all.adjusted,rep("cp.all.adjusted",length(cp.all.adjusted))),
# cbind(cp.all.unadjusted,rep("cp.all.unadjusted",length(cp.all.unadjusted))))
# 
# colnames(yy)<-c("y","grp")
# 
# yy<-as.data.frame(yy)
# head(yy)
# boxplot(as.numeric(as.character(y))~grp,data=yy)
# 
# Output_file="/Volumes/Bioinformatics$/Aimin_project/ToGaoZhen/In_ad_not_un.xls"
# WriteGoToTable(Example.Go.adjusted.by.exon$GO.selected[index1,],Output_file)
# 
# Output_file="/Volumes/Bioinformatics$/Aimin_project/ToGaoZhen/In_un_not_ad.xls"
# WriteGoToTable(Example.Go.unadjusted$GO.selected[index2,],Output_file)







