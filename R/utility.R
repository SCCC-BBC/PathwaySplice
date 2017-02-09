# Utility functions for PathwaySplice

list_to_df <- function(list_for_df) {
    list_for_df <- as.list(list_for_df)
    
    nm <- names(list_for_df)
    if (is.null(nm)) 
        nm <- seq_along(list_for_df)
    
    df <- data.frame(name = nm, stringsAsFactors = FALSE)
    df$value <- unname(list_for_df)
    df
}

reversemapping = function(map) {
    tmp = unlist(map, use.names = FALSE)
    names(tmp) = rep(names(map), times = as.numeric(summary(map)[, 
        1]))
    return(split(names(tmp), as.vector(tmp)))
}

ReformatData <- function(re.PJ.gene.based) {
    re <- pData(re.PJ.gene.based)
    
    no.re.testable.index <- which(as.character(re$mostSigID) == 
        "character(0)")
    re2 <- re[-no.re.testable.index, ]
    
    All.gene.id.based.on.sub_feature <- unique(unlist(strsplit(re2[, 
        1], "\\+")))
    All.gene.id.index <- rep(0, length(All.gene.id.based.on.sub_feature))
    names(All.gene.id.index) = All.gene.id.based.on.sub_feature
    
    
    reformat.gene.p <- do.call(rbind, sapply(All.gene.id.based.on.sub_feature, 
        function(u, re2) {
            x <- re2[grep(u, re2[, 1]), -1]
            x <- as.data.frame(t(x))
            # colnames(x)<-colnames(Data4Goterm) x
        }, re2))
    
    re3 <- as.data.frame(reformat.gene.p)
    re3 <- cbind(All.gene.id.based.on.sub_feature, re3)
    colnames(re3)[1] = "geneID"
    
    return(re3)
    
}

# @examples dir.name='/media/H_driver/2016/Yang/MACS/MACS/'
# reformatPath(dir.name)

reformatPath <- function(dir.name) {
    CheckOPS = Sys.info()[["sysname"]]
    
    if (CheckOPS == "Darwin") {
        temp = unlist(strsplit(dir.name, split = "\\/"))
        
        #unlist(strsplit('/media/H_driver/Aimin_project/',split="\\/"))
        
        if(temp[3]=="H_driver"){
        temp[2] = "Volumes"
        temp[3] = "Bioinformatics$"
        dir.name = paste0(paste0(temp, collapse = "/"), "/")
       }
    }
    
    return(dir.name)
}


heatmap_wPCA = function(Data,g_level = NULL) {
  
  Data.pca = prcomp(t(Data))
  hmcol<-rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
  
  if(is.null(g_level)) {
    type_level = 1:ncol(Data)
    col_level = "black"
    
    with(data.frame(Data.pca$x), scatter3D(PC1, PC2, PC3, colvar = NULL, type="h",
                                           ticktype = "detailed", bty="b2", cex=1,
                                           xlab="PC 1",	ylab="PC 2",zlab="PC 3", theta = 40, phi = 40, pch=type_level,
                                           col=col_level,
                                           main = "Principal component analysis"))
    
    
    
    with(data.frame(Data.pca$x), text3D(x=PC1, y=PC2,
                                        z=PC3, colnames(Data), col = "black", add=TRUE, colkey = FALSE, cex=0.5))
    
    
    
    
  } else {
    type_level = 1:ncol(Data)
    TEMP = factor(g_level)
    uniq_label =  levels(TEMP)
    levels(TEMP) = hmcol[ceiling(seq(length.out=length(levels(TEMP)),from=1,to=256))]
    col_level = as.character(TEMP)
    uniq_col = levels(TEMP)
    
    Data.pca = prcomp(t(Data))
    with(data.frame(Data.pca$x), scatter3D(PC1, PC2, PC3, colvar = NULL, type="h",
                                           ticktype = "detailed", bty="b2", cex=1,
                                           xlab="PC 1",	ylab="PC 2",zlab="PC 3", theta = 40, phi = 40, pch=type_level,
                                           col=col_level,
                                           main = "Principal component analysis"))
    
    legend("topright", legend = uniq_label, pch=type_level,
           col = uniq_col,
           cex=1, inset=c(0.02))
    
    with(data.frame(Data.pca$x), text3D(x=PC1, y=PC2,
                                        z=PC3, colnames(Data), col = "black", add=TRUE, colkey = FALSE, cex=0.5))
  }
}

WriteGoToTable <- function(GO_re,Output_file) {
  dataset2<- GO_re
  dataset2[sapply(dataset2, is.list)] <-
    sapply(dataset2[sapply(dataset2, is.list)],
           function(x)sapply(x, function(y) paste(unlist(y),collapse=", ") ) )
  
  write.table(dataset2,file=Output_file,row.names = FALSE,quote=FALSE,sep="\t")
}


#' PostProcessGO
#'
#' @param n.go 
#' @param adjusted 
#' @param unadjuasted 
#' @param venn.dir 
#' @param boxplot.dir 
#' @param In.ad.not.un.file 
#' @param In.un.not.ad.file 
#'
#' @return null
#' 
#'
#' @examples
#' 
#' PostProcessGO(25,Example.cp.adjusted.by.exon,Example.cp.unadjusted,
#' "/Volumes/Bioinformatics$/Aimin_project/ToGaoZhen/","/Volumes/Bioinformatics$/Aimin_project/ToGaoZhen/",
#' In_ad_not_un.xls","In_un_not_ad.xls")
#' 
#' @export
#' 
#' 
#' 
#' 
PostProcessGO <- function(n.go,adjusted,unadjuasted,venn.dir,boxplot.dir,In.ad.not.un.file,In.un.not.ad.file) {
  n=n.go
  
  Example.Go.adjusted.by.exon<-adjusted
  Example.Go.unadjusted<-unadjuasted
  
  adjusted<-Example.Go.adjusted.by.exon$GO.selected[1:n,1]
  unadjusted<-Example.Go.unadjusted$GO.selected[1:n,1]
  
  re<-list(adjusted=adjusted,unadjusted=unadjusted)
  
  venn.plot <- venn.diagram(
    x = re[c(1,2)],
    filename = paste0(venn.dir,"/",names(re)[1],"_",names(re)[2],"_overlap_venn.tiff"),
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
  In.ad.not.un<-Example.Go.adjusted.by.exon$GO.selected[index1,]$Ave_value_all_gene
  
  index2<-match(In.unadjusted.not.in.adjusted,Example.Go.unadjusted$GO.selected$category)
  In.un.not.ad<-Example.Go.unadjusted$GO.selected[index2,]$Ave_value_all_gene
  
  xx<-cbind(unlist(In.ad.not.un),unlist(In.un.not.ad))
  
  colnames(xx)=c("In.ad.not.un","In.un.not.ad")
  
  #boxplot(xx)
  
  #cbind(Example.Go.adjusted.by.exon$GO.selected[index1,1],Example.Go.unadjusted$GO.selected[index2,1])
  
  #In.ad.not.un<-xx[,1]
  #In.un.not.ad<-xx[,2]
  
  cp.top.adjusted.25<-unlist(Example.Go.adjusted.by.exon$GO.selected[1:n,]$Ave_value_all_gene)
  cp.top.unadjusted.25<-unlist(Example.Go.unadjusted$GO.selected[1:n,]$Ave_value_all_gene)
  
  cp.all.adjusted<-unlist(Example.Go.adjusted.by.exon$GO.selected$Ave_value_all_gene)
  cp.all.unadjusted<-unlist(Example.Go.unadjusted$GO.selected$Ave_value_all_gene)
  
  yy<-rbind(cbind(xx[,1],rep("In.ad.not.un",length(xx[,1]))),
            cbind(xx[,2],rep("In.un.not.ad",length(xx[,2]))),
            cbind(cp.top.adjusted.25,rep("cp.top.adjusted.25",length(cp.top.adjusted.25))),
            cbind(cp.top.unadjusted.25,rep("cp.top.unadjusted.25",length(cp.top.unadjusted.25))),
            cbind(cp.all.adjusted,rep("cp.all.adjusted",length(cp.all.adjusted))),
            cbind(cp.all.unadjusted,rep("cp.all.unadjusted",length(cp.all.unadjusted))))
  
  colnames(yy)<-c("y","grp")
  
  yy<-as.data.frame(yy)
  #head(yy)
  
  png(paste0(boxplot.dir,"/","boxplot.png"))
  boxplot(as.numeric(as.character(y))~grp,data=yy)
  dev.off()
    
  Output_file=paste0(boxplot.dir,"/",In.ad.not.un.file)
  WriteGoToTable(Example.Go.adjusted.by.exon$GO.selected[index1,],Output_file)
  
  Output_file=paste0(boxplot.dir,"/",In.un.not.ad.file)
  WriteGoToTable(Example.Go.unadjusted$GO.selected[index2,],Output_file)
}

















