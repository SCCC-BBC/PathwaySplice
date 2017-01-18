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

