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

reversemapping <- function(map) {
    tmp <- unlist(map, use.names = FALSE)
    names(tmp) <- rep(names(map), times = as.numeric(summary(map)[, 1]))
    return(split(names(tmp), as.vector(tmp)))
}

reformatdata <- function(re.PJ.gene.based) {
    re <- pData(re.PJ.gene.based)
    
    no.re.testable.index <- which(as.character(re$mostSigID) == "character(0)")
    re2 <- re[-no.re.testable.index, ]
    
    All.gene.id.based.on.sub_feature <- unique(unlist(strsplit(re2[, 1], "\\+")))
    All.gene.id.index <- rep(0, length(All.gene.id.based.on.sub_feature))
    names(All.gene.id.index) <- All.gene.id.based.on.sub_feature
    
    
    reformat.gene.p <- do.call(rbind, sapply(All.gene.id.based.on.sub_feature, function(u, 
        re2) {
        x <- re2[grep(u, re2[, 1]), -1]
        x <- as.data.frame(t(x))
        # colnames(x)<-colnames(Data4Goterm) x
    }, re2))
    
    re3 <- as.data.frame(reformat.gene.p)
    re3 <- cbind(All.gene.id.based.on.sub_feature, re3)
    colnames(re3)[1] <- "geneID"
    
    return(re3)
    
}

#' reformatpath
#' 
#' @param dir.name Directory name to be converted 
#' @return A converted directory
#' 
#' @examples
#' dir.name <- '/media/H_driver/2016/Yang/MACS/MACS/'
#' reformatpath(dir.name)
#' 
#' @export 

reformatpath <- function(dir.name) {
    CheckOPS <- Sys.info()[["sysname"]]
    
    if (CheckOPS == "Darwin") {
        temp <- unlist(strsplit(dir.name, split = "\\/"))
        
        if (!is.na(temp[3] == "H_driver")) {
            if (temp[3] == "H_driver") {
                temp[2] <- "Volumes"
                temp[3] <- "Bioinformatics$"
                dir.name <- file.path(paste0(temp, collapse = "/"))
            }
        }
        
    }
    
    return(dir.name)
}

heatmap_wPCA = function(Data, g_level = NULL) {
    
    Data.pca = prcomp(t(Data))
    hmcol <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
    
    if (is.null(g_level)) {
        type_level = 1:ncol(Data)
        col_level = "black"
        
        with(data.frame(Data.pca$x), scatter3D(PC1, PC2, PC3, colvar = NULL, type = "h", 
            ticktype = "detailed", bty = "b2", cex = 1, xlab = "PC 1", ylab = "PC 2", 
            zlab = "PC 3", theta = 40, phi = 40, pch = type_level, col = col_level, 
            main = "Principal component analysis"))
        
        
        
        with(data.frame(Data.pca$x), text3D(x = PC1, y = PC2, z = PC3, colnames(Data), 
            col = "black", add = TRUE, colkey = FALSE, cex = 0.5))
        
        
        
        
    } else {
        type_level = 1:ncol(Data)
        TEMP = factor(g_level)
        uniq_label = levels(TEMP)
        levels(TEMP) = hmcol[ceiling(seq(length.out = length(levels(TEMP)), from = 1, 
            to = 256))]
        col_level = as.character(TEMP)
        uniq_col = levels(TEMP)
        
        Data.pca = prcomp(t(Data))
        with(data.frame(Data.pca$x), scatter3D(PC1, PC2, PC3, colvar = NULL, type = "h", 
            ticktype = "detailed", bty = "b2", cex = 1, xlab = "PC 1", ylab = "PC 2", 
            zlab = "PC 3", theta = 40, phi = 40, pch = type_level, col = col_level, 
            main = "Principal component analysis"))
        
        legend("topright", legend = uniq_label, pch = type_level, col = uniq_col, 
            cex = 1, inset = c(0.02))
        
        with(data.frame(Data.pca$x), text3D(x = PC1, y = PC2, z = PC3, colnames(Data), 
            col = "black", add = TRUE, colkey = FALSE, cex = 0.5))
    }
}

writegototable <- function(GO_re, Output_file) {
    dataset2 <- GO_re
    dataset2[sapply(dataset2, is.list)] <- sapply(dataset2[sapply(dataset2, is.list)], 
        function(x) sapply(x, function(y) paste(unlist(y), collapse = ", ")))
    
    write.table(dataset2, file = Output_file, row.names = FALSE, quote = FALSE, sep = "\t")
}

#' postprocessgo
#'
#' @param n.go Number of gene sets
#' @param adjusted Adjusted result 
#' @param unadjuasted Unadjusted result 
#' @param venn.dir Path for outputing venn 
#' @param boxplot.dir Path for outputing boxplot
#' @param type.boxplot Get boxplot for 5 categories or 3 categories
#'        6 categories: 'All.adjusted','All.unadjusted',
#'                      'Top25.adjusted','Top25.unadjusted',
#'                      'In_ad_not_un','In_un_not_ad'
#'        3 categories: 'All','Top25.adjusted','Top25.unadjusted
#' @param In.ad.not.un.file File name for outputing adjused but not
#'        in unadjusted when using the selected gene sets     
#' @param In.un.not.ad.file File name for outputing unadjused but not
#'        in adjusted when using the selected gene sets  
#'
#' @return null
#' 
#' @examples
#' 
#' dir.name <- system.file('extdata', package='PathwaySplice')
#' canonical.pathway.file <- '10.cp.gmt.txt'
#' res <- gmtgene2cat(dir.name,canonical.pathway.file,'local',genomeID='hg19')

#' res1 <- runpathwaysplice(tiny.data,adjust='exon_SJ',sub_feature='E',
#' 0.05,genomeID='hg19',geneID='ensGene',gene2cat=res,method='Wallenius')

#' res2 <- runpathwaysplice(tiny.data,adjust='exon_SJ',sub_feature='E',
#' 0.05,genomeID='hg19',geneID='ensGene',gene2cat=res,method='Hypergeometric')

#' dir.name <- tempdir()
#' output.dir <- file.path(dir.name,"OutputPostAnalysis")
#' 
#' output.file.name.1 <- "In_ad_not_un.xls"
#' output.file.name.2 <- "In_un_not_ad.xls"
#' res3 <- postprocessgo(4,res1,res2,output.dir,output.dir,
#' type.boxplot='Only3',output.file.name.1,output.file.name.2)

#' @export

postprocessgo <- function(n.go, adjusted, unadjuasted, venn.dir, boxplot.dir, type.boxplot = c("All", 
    "Only3"), In.ad.not.un.file, In.un.not.ad.file) {
    
    if (!dir.exists(venn.dir)) {
        dir.create(venn.dir)
    }
    
    if (!dir.exists(boxplot.dir)) {
        dir.create(boxplot.dir)
    }
    
    n <- n.go
    
    Example.Go.adjusted.by.exon <- adjusted
    Example.Go.unadjusted <- unadjuasted
    
    
    if (dim(Example.Go.adjusted.by.exon$GO.selected)[1]>=n&&
        dim(Example.Go.unadjusted$GO.selected)[1]>=n){
      
    adjusted <- Example.Go.adjusted.by.exon$GO.selected[1:n, 1]
    unadjusted <- Example.Go.unadjusted$GO.selected[1:n, 1]
    
    re <- list(adjusted = adjusted, unadjusted = unadjusted)
    
    venn.plot <- venn.diagram(x = re[c(1, 2)], filename = file.path(venn.dir, paste0(names(re)[1], 
        "_", names(re)[2], "_overlap_venn.tiff")), height = 3000, width = 3500, resolution = 1000, 
        col = "black", lty = "dotted", lwd = 1, fill = c("red", "blue"), alpha = 0.5, 
        label.col = c(rep("black", 3)), cex = 0.5, fontfamily = "serif", fontface = "bold", 
        cat.col = c("red", "blue"), cat.cex = 0.5, cat.pos = 0.5, cat.dist = 0.05, 
        cat.fontfamily = "serif")
    
    # boxplot
    common <- intersect(unadjusted, adjusted)
    
    In.unadjusted.not.in.adjusted <- setdiff(unadjusted, common)
    In.adjusted.not.in.unadjusted <- setdiff(adjusted, common)
    
    if (length(In.unadjusted.not.in.adjusted) != 0 && length(In.adjusted.not.in.unadjusted) != 
        0) {
        index1 <- match(In.adjusted.not.in.unadjusted, Example.Go.adjusted.by.exon$GO.selected$category)
        In.ad.not.un <- Example.Go.adjusted.by.exon$GO.selected[index1, ]$Ave_value_all_gene
        
        yy <- cbind(Example.Go.unadjusted$GO.selected[index1, ]$rank.value.by.over_represented_pvalue, 
            Example.Go.adjusted.by.exon$GO.selected[index1, ]$rank.value.by.over_represented_pvalue)
        
        
        index2 <- match(In.unadjusted.not.in.adjusted, Example.Go.unadjusted$GO.selected$category)
        In.un.not.ad <- Example.Go.unadjusted$GO.selected[index2, ]$Ave_value_all_gene
        
        yyy <- cbind(Example.Go.unadjusted$GO.selected[index2, ]$rank.value.by.over_represented_pvalue, 
            Example.Go.adjusted.by.exon$GO.selected[index2, ]$rank.value.by.over_represented_pvalue)
        
        rre <- list(yy = yy, yyy = yyy)
        
        xx <- cbind(unlist(In.ad.not.un), unlist(In.un.not.ad))
        
        colnames(xx) <- c("In.ad.not.un", "In.un.not.ad")
        
        cp.top.adjusted.25 <- unlist(Example.Go.adjusted.by.exon$GO.selected[1:n, 
            ]$Ave_value_all_gene)
        cp.top.unadjusted.25 <- unlist(Example.Go.unadjusted$GO.selected[1:n, ]$Ave_value_all_gene)
        
        cp.all.adjusted <- unlist(Example.Go.adjusted.by.exon$GO.selected$Ave_value_all_gene)
        cp.all.unadjusted <- unlist(Example.Go.unadjusted$GO.selected$Ave_value_all_gene)
        
        type.boxplot <- match.arg(type.boxplot)
        
        switch(type.boxplot, Only3 = {
            yy <- rbind(cbind(cp.top.adjusted.25, rep("Adjusted_25", length(cp.top.adjusted.25))), 
                cbind(cp.top.unadjusted.25, rep("Unadjusted_25", length(cp.top.unadjusted.25))), 
                cbind(cp.all.unadjusted, rep("All", length(cp.all.unadjusted))))
            colnames(yy) <- c("y", "grp")
            yy <- as.data.frame(yy)
            yy$grp <- factor(yy$grp)
            yy$grp <- factor(yy$grp, levels = levels(yy$grp)[c(2, 1, 3)])
            png(file.path(boxplot.dir, "boxplot.png"))
            boxplot(as.numeric(as.character(y)) ~ grp, data = yy)
            dev.off()
        }, {
            yy <- rbind(cbind(xx[, 1], rep("In.ad.not.un", length(xx[, 1]))), cbind(xx[, 
                2], rep("In.un.not.ad", length(xx[, 2]))), cbind(cp.top.adjusted.25, 
                rep("cp.top.adjusted.25", length(cp.top.adjusted.25))), cbind(cp.top.unadjusted.25, 
                rep("cp.top.unadjusted.25", length(cp.top.unadjusted.25))), cbind(cp.all.adjusted, 
                rep("cp.all.adjusted", length(cp.all.adjusted))), cbind(cp.all.unadjusted, 
                rep("cp.all.unadjusted", length(cp.all.unadjusted))))
            colnames(yy) <- c("y", "grp")
            yy <- as.data.frame(yy)
            png(file.path(boxplot.dir, "boxplot.png"))
            boxplot(as.numeric(as.character(y)) ~ grp, data = yy)
            dev.off()
        })
        
        Output_file <- file.path(boxplot.dir, In.ad.not.un.file)
        writegototable(Example.Go.adjusted.by.exon$GO.selected[index1, ], Output_file)
        
        Output_file <- file.path(boxplot.dir, In.un.not.ad.file)
        writegototable(Example.Go.unadjusted$GO.selected[index2, ], Output_file)
        
        return(rre)
        
    } else {
        
        if (length(In.unadjusted.not.in.adjusted) == 0) {
            cat("there is no gene sets in unadjusted resutls but not in adjusted resutls\n")
        }
        
        cat("\n")
        
        if (length(In.adjusted.not.in.unadjusted) == 0) {
            cat("there is no gene sets in adjusted resutls but not in unadjusted resutls\n")
        }
        
        cat("\n")
        
    }
    } else {
    
      cat("The enriched gene sets is less than", n,"\n")
      
    }

}

match2Genome <- function(genome_id) {
    
    ah <- AnnotationHub()
    switch(genome_id, hg38 = {
        hs <- query(ah, c("Ensembl", "GRCh38", "Homo sapiens"))
        res <- hs[["AH53211"]]
        res <- genes(res, columns = c("gene_name"))
        xxx <- mcols(res)
        yyy <- xxx
    }, hg19 = {
        edb <- org.Hs.eg.db
        entrezid <- keys(edb, keytype = "ENTREZID")
        suppressMessages(xxx <- select(edb, keys = entrezid, columns = c("ENSEMBL", 
            "SYMBOL")))
        yyy <- xxx[, c(3, 2, 1)]
    }, {
        hs <- query(ah, c("Ensembl", "GRCm38", "Mus Musculus"))
        res <- hs[["AH53222"]]
        res <- genes(res, columns = c("gene_name"))
        xxx <- S4Vectors::mcols(res)
        yyy <- xxx
    })
    return(yyy)
}
