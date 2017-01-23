library(VennDiagram)

#Venn
n=25
adjusted<-Example.Go.adjusted.by.exon$GO.selected[1:n,1]
unadjusted<-Example.Go.unadjusted$GO.selected[1:n,1]

re<-list(adjusted=adjusted,unadjusted=unadjusted)

venn.plot <- venn.diagram(
  x = re[c(1,2)],
  filename = paste0(getwd(),"/",names(re)[1],"_",names(re)[2],"_overlap_venn.tiff"),
  #filename=NULL,
  height = 3000,
  width = 3500,
  resolution = 1000,
  col = "black",
  lty = "dotted",
  lwd = 1,
  fill = c("red","blue"),
  alpha = 0.50,
  label.col = c(rep("black",3)),
  cex = 0.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("red","blue"),
  cat.cex = 0.5,
  cat.pos = 0.5,
  cat.dist = 0.05,
  cat.fontfamily = "serif"
)


#boxplot

common<-intersect(unadjusted,adjusted)

In.unadjusted.not.in.adjusted<-setdiff(unadjusted,common)
In.adjusted.not.in.unadjusted<-setdiff(adjusted,common)


length(In.unadjusted.not.in.adjusted)
length(In.adjusted.not.in.unadjusted)
length(common)

index1<-match(In.adjusted.not.in.unadjusted,Example.Go.adjusted.by.exon$GO.selected$category)
In.ad.not.un<-Example.Go.adjusted.by.exon$GO.selected[index,]$Ave_value_all_gene

index2<-match(In.unadjusted.not.in.adjusted,Example.Go.unadjusted$GO.selected$category)

In.un.not.ad<-Example.Go.unadjusted$GO.selected[index2,]$Ave_value_all_gene

xx<-cbind(unlist(In.ad.not.un),unlist(In.un.not.ad))

colnames(xx)=c("In.ad.not.un","In.un.not.ad")

boxplot(xx)

