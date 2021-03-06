# To generate this file:
# 
# library(formatR)
# tidy_source(source="R/Run_pathwaysplice.R", comment = FALSE,file="inst/bin/Simple.r")
# 
makeGeneTable <- function(feature.table, sig.threshold = 0.05, stat = "pvalue") {
    min.pval <- aggregate(pvalue ~ geneID, data = feature.table, FUN = min)
    n.feature <- as.data.frame(table(feature.table$geneID))
    both <- merge(x = min.pval, y = n.feature, by.x = "geneID", by.y = "Var1")
    both$fdr <- p.adjust(both$pvalue, method = "fdr")
    if (stat == "pvalue") {
        both$sig.gene <- ifelse(both$pvalue < sig.threshold, 1, 0)
    }
    if (stat == "fdr") {
        both$sig.gene <- ifelse(both$fdr < sig.threshold, 1, 0)
    }
    names(both) <- sub("pvalue", "geneWisePvalue", names(both))
    names(both) <- sub("Freq", "numFeature", names(both))
    return(both)
}
lrTestBias <- function(genewise.table, boxplot.width = 0.1) {
    mydata <- genewise.table
    DE <- ifelse(mydata$sig.gene == 1, 1, 0)
    mydata.2 <- cbind(mydata, DE)
    if (var(mydata.2$numFeature) != 0) {
        mylogit.2 <- glm(DE ~ as.numeric(numFeature), data = mydata.2, family = "binomial")
        re <- summary(mylogit.2)
        pvalue <- re$coefficients[2, 4]
        p.value <- format.pval(pvalue, eps = 1e-04, digits = 2)
        boxplot(mydata.2$numFeature ~ mydata.2$DE, boxwex = boxplot.width, ylab = "Number of features", col = "lightgray", ylim = c(min(mydata.2$numFeature), 
            max(mydata.2$numFeature)), names = c("non-significant genes", "significant genes"))
        text(x = 2, y = max(mydata.2$numFeature) - 1, labels = c("", paste0("P-value from logistic regression ", p.value)))
        
    }
    else {
        cat("There are no variations on the number of features\n")
    }
}
runPathwaySplice <- function(genewise.table, genome, id, gene2cat = NULL, test.cats = c("GO:CC", "GO:BP", "GO:MF"), go.size.limit = c(10, 
    200), method = "Wallenius", repcnt = 2000, use.genes.without.cat = FALSE, binsize = "auto", output.file = tempfile()) {
    x <- genewise.table$sig.gene
    names(x) <- genewise.table$geneID
    pwf <- nullpSplice(x, genome, id, bias.data = genewise.table$numFeature, plot.fit = TRUE, binsize)
    CatDE <- pathwaysplice(pwf, genome = genome, id = id, gene2cat = gene2cat, test.cats = test.cats, go.size.limit = go.size.limit, method = method, 
        repcnt = repcnt, use.genes.without.cat = use.genes.without.cat)
    res1 <- getStaisitcs4Go(CatDE, genewise.table)
    res2 <- reformatPathwayOut(res1)
    res3 <- within(res2, rm("Ave_value_DE"))
    writeTibble(res3, output.file)
    res3
}
enrichmentMap <- function(pathway.res, n = 50, fixed = TRUE, node.label.font = 1, similarity.threshold, scaling.factor = 1, output.file.dir = tempdir(), 
    label.node.by.index = FALSE, ...) {
    if (!dir.exists(output.file.dir)) {
        dir.create(output.file.dir, recursive = TRUE)
    }
    pathway.res <- as.data.frame(pathway.res)
    GO.name <- pathway.res$gene_set
    temp <- pathway.res$SIGgene_ensembl
    names(temp) <- GO.name
    x <- pathway.res
    geneSets <- temp
    y <- as.data.frame(x)
    if (any(grep("^GO:", y$gene_set))) {
        nodename <- paste0(y$description, ":", y$numSIGInCat)
        if (label.node.by.index == TRUE) {
            nodename.index <- seq(1, length(nodename))
            output.text <- as.data.frame(cbind(nodename.index, nodename))[1:n, ]
            nodename <- nodename.index
            colnames(output.text) <- c("index", "name")
            write.table(output.text, file = file.path(output.file.dir, "enrichmap_GO.xls"), quote = FALSE, col.names = TRUE, row.names = FALSE, 
                sep = "\t")
        }
    }
    else {
        nodename <- paste0(y$gene_set, ":", y$numSIGInCat)
        if (label.node.by.index == TRUE) {
            nodename.index <- seq(1, length(nodename))
            output.text <- as.data.frame(cbind(nodename.index, nodename))[1:n, ]
            nodename <- nodename.index
            colnames(output.text) <- c("index", "name")
            write.table(output.text, file = file.path(output.file.dir, "enrichmap_pathway.xls"), quote = FALSE, col.names = TRUE, row.names = FALSE, 
                sep = "\t")
        }
    }
    if (nrow(y) < n) {
        n <- nrow(y)
    }
    y <- y[1:n, ]
    if (n == 0) {
        stop("no enriched term found...")
    }
    else if (n == 1) {
        g <- graph.empty(0, directed = FALSE)
        g <- add_vertices(g, nv = 1)
        igraph::V(g)$name <- nodename
        igraph::V(g)$color <- "red"
    }
    else {
        pvalue <- as.numeric(y$over_represented_pvalue)
        id <- y[, 1]
        geneSets <- geneSets[id]
        n <- nrow(y)
        w <- matrix(NA, nrow = n, ncol = n)
        colnames(w) <- rownames(w) <- nodename[1:n]
        for (i in 1:n) {
            for (j in i:n) {
                w[i, j] <- overlap_ratio(geneSets[id[i]], geneSets[id[j]])
            }
        }
        wd <- melt(w)
        wd <- wd[wd[, 1] != wd[, 2], ]
        wd <- wd[!is.na(wd[, 3]), ]
        g <- graph.data.frame(wd[, -3], directed = FALSE)
        igraph::E(g)$width <- sqrt(wd[, 3] * 5) * scaling.factor
        g <- delete.edges(g, igraph::E(g)[wd[, 3] < similarity.threshold])
        idx <- unlist(sapply(igraph::V(g)$name, function(x) which(x == nodename[1:n])))
        cols <- color_scale("red", "#E5C494")
        igraph::V(g)$color <- cols[sapply(pvalue, getIdx, min(pvalue), max(pvalue))]
        Edata <- as.data.frame(get.edgelist(g))
        Edata$edgewidth <- igraph::E(g)$width
        Vdata <- data.frame(pathway = igraph::V(g)$name, color = igraph::V(g)$color)
        map_data <- list(edge_data = Edata, node_data = Vdata)
        cnt <- as.integer(y$numSIGInCat)
        names(cnt) <- nodename[1:n]
        cnt2 <- cnt[igraph::V(g)$name]
        igraph::V(g)$size <- cnt2/sum(cnt2) * 100
    }
    netplot(g, node.label.font = node.label.font, node.label.color = "black", fixed = fixed, ...)
    write.graph(g, file.path(output.file.dir, "network.layout.for.cytoscape.gml"), format = "gml")
    invisible(g)
    re2 <- map_data
    return(re2)
}
gmtGene2Cat <- function(pathway.file, gene.anno.file = NULL, genomeID = c("mm10", "hg19", "hg38")) {
    gmt_input_file <- pathway.file
    gene.2.cat.gmt <- gene2cat2(gmt_input_file)
    names.gene.gmt <- as.data.frame(names(gene.2.cat.gmt))
    colnames(names.gene.gmt) <- "gene_id"
    if (!is.null(gene.anno.file)) {
        gene_anno_file <- gene.anno.file
        gene.id.conversion <- read.csv(gene_anno_file)
    }
    else {
        gene.id.conversion <- match.arg(genomeID)
    }
    xxx <- match2Genome(gene.id.conversion)
    names.gene.gmt.2 <- match(names.gene.gmt$gene_id, xxx[, 1])
    gene.id.conversion.2 <- xxx[names.gene.gmt.2, ]
    gene.2.cat.gmt.2 <- gene.2.cat.gmt
    names(gene.2.cat.gmt.2) <- gene.id.conversion.2[, 2]
    gene.2.cat.gmt.2
}
.onAttach <- function(libname, pkgname) {
    if (.Platform$OS.type == "windows" && .Platform$GUI == "Rgui") {
        winMenuAddItem("Vignettes", "PathwaySplice", "shell.exec(system.file(\"doc\",\"PathwaySplice.pdf\",package=\"PathwaySplice\"))")
    }
}
.ORG_PACKAGES = paste("org.", c("Ag.eg", "At.tair", "Bt.eg", "Ce.eg", "Cf.eg", "Dm.eg", "Dr.eg", "EcK12.eg", "EcSakai.eg", "Gg.eg", "Hs.eg", 
    "Mm.eg", "Mmu.eg", "Pf.plasmo", "Pt.eg", "Rn.eg", "Sc.sgd", "Ss.eg", "Xl.eg"), sep = "")
names(.ORG_PACKAGES) = c("anoGam", "Arabidopsis", "bosTau", "ce", "canFam", "dm", "danRer", "E. coli K12", "E. coli Sakai", "galGal", "hg", 
    "mm", "rheMac", "Malaria", "panTro", "rn", "sacCer", "susScr", "xenTro")
.ID_MAP = c("eg", "eg", "ENSEMBL", "SYMBOL", "sgd", "plasmo", "tair")
names(.ID_MAP) = c("knownGene", "refGene", "ensGene", "geneSymbol", "sgd", "plasmo", "tair")
.ORG_GOMAP_FUNCTION = c("GO2ALLEGS", "GO2ALLTAIRS", "GO2ALLORFS", "GO2ALLORFS")
names(.ORG_GOMAP_FUNCTION) = c("default", "org.At.tair", "org.Pf.plasmo", "org.Sc.sgd")
.TXDB_ORGS = c("ce6", "dm3", "hg18", "hg19", "hg38", "mm10", "mm9", "rn4", "rn5", "sacCer2", "sacCer3")
compareResults <- function(n.go, adjusted, unadjusted, gene.based.table, output.dir = tempdir(), type.boxplot = c("All", "Only3")) {
    if (!dir.exists(output.dir)) {
        dir.create(output.dir, recursive = TRUE)
    }
    n <- n.go
    example.go.adjusted.by.exon <- as.data.frame(adjusted)
    example.go.unadjusted <- as.data.frame(unadjusted)
    if (is.na(example.go.adjusted.by.exon) || is.na(example.go.unadjusted)) {
        cat("One of results is empty\n\n")
        return()
    }
    else {
        if (dim(example.go.adjusted.by.exon)[1] >= n && dim(example.go.unadjusted)[1] >= n) {
            adjusted <- example.go.adjusted.by.exon[1:n, 1]
            unadjusted <- example.go.unadjusted[1:n, 1]
            re <- list(adjusted = adjusted, unadjusted = unadjusted)
            vp <- venn.diagram(re, fill = c("red", "blue"), cat.col = c("red", "blue"), alpha = 0.3, filename = file.path(output.dir, paste0(names(re)[1], 
                "_", names(re)[2], "_overlap_venn.tiff")))
            common <- intersect(unadjusted, adjusted)
            In.unadjusted.not.in.adjusted <- setdiff(unadjusted, common)
            In.adjusted.not.in.unadjusted <- setdiff(adjusted, common)
            if (length(In.unadjusted.not.in.adjusted) != 0 && length(In.adjusted.not.in.unadjusted) != 0) {
                index1 <- match(In.adjusted.not.in.unadjusted, example.go.adjusted.by.exon$gene_set)
                index1.name <- unique(unlist(example.go.adjusted.by.exon[index1, ]$All_gene_ensembl))
                In.ad.not.un <- gene.based.table[match(index1.name, gene.based.table$geneID), ]$numFeature
                yy <- cbind(example.go.unadjusted[index1, ]$rank.value.by.over_represented_pvalue, example.go.adjusted.by.exon[index1, 
                  ]$rank.value.by.over_represented_pvalue)
                index2 <- match(In.unadjusted.not.in.adjusted, example.go.unadjusted$gene_set)
                index2.name <- unique(unlist(example.go.unadjusted[index2, ]$All_gene_ensembl))
                In.un.not.ad <- gene.based.table[match(index2.name, gene.based.table$geneID), ]$numFeature
                yyy <- cbind(example.go.unadjusted[index2, ]$rank.value.by.over_represented_pvalue, example.go.adjusted.by.exon[index2, 
                  ]$rank.value.by.over_represented_pvalue)
                cp.topN.adjusted.name <- unique(unlist(example.go.adjusted.by.exon[1:n, ]$All_gene_ensembl))
                cp.topN.adjusted <- gene.based.table[match(cp.topN.adjusted.name, gene.based.table$geneID), ]$numFeature
                cp.topN.unadjusted.name <- unique(unlist(example.go.unadjusted[1:n, ]$All_gene_ensembl))
                cp.topN.unadjusted <- gene.based.table[match(cp.topN.unadjusted.name, gene.based.table$geneID), ]$numFeature
                adjusted.name <- unique(unlist(example.go.adjusted.by.exon$All_gene_ensembl))
                cp.all.adjusted <- gene.based.table[match(adjusted.name, gene.based.table$geneID), ]$numFeature
                unadjusted.name <- unique(unlist(example.go.unadjusted$All_gene_ensembl))
                cp.all.unadjusted <- gene.based.table[match(unadjusted.name, gene.based.table$geneID), ]$numFeature
                type.boxplot <- match.arg(type.boxplot)
                switch(type.boxplot, Only3 = {
                  yy <- rbind(cbind(cp.topN.adjusted, rep(paste0("Adjusted_", n.go), length(cp.topN.adjusted))), cbind(cp.topN.unadjusted, 
                    rep(paste0("Unadjusted_", n.go), length(cp.topN.unadjusted))), cbind(cp.all.unadjusted, rep("All", length(cp.all.unadjusted))))
                  colnames(yy) <- c("y", "grp")
                  yy <- as.data.frame(yy)
                  yy$grp <- factor(yy$grp)
                  yy$grp <- factor(yy$grp, levels = levels(yy$grp)[c(2, 1, 3)])
                  yy$y <- as.numeric(as.character(yy$y))
                  colnames(yy) <- c("numFeature", "category")
                  m = list(l = 200, r = 5, b = 5, t = 5, pad = 4)
                  p <- plot_ly(yy, x = ~numFeature, color = ~category, type = "box") %>% layout(margin = m)
                  htmlwidgets::saveWidget(p, file.path(output.dir, "boxplot.html"))
                }, {
                  yy <- rbind(cbind(In.ad.not.un, rep("adjusted.only", length(In.ad.not.un))), cbind(In.un.not.ad, rep("unadjusted.only", 
                    length(In.un.not.ad))), cbind(cp.topN.adjusted, rep(paste0("top.adjusted.", n.go), length(cp.topN.adjusted))), cbind(cp.topN.unadjusted, 
                    rep(paste0("top.unadjusted.", n.go), length(cp.topN.unadjusted))), cbind(cp.all.adjusted, rep("all.adjusted", length(cp.all.adjusted))), 
                    cbind(cp.all.unadjusted, rep("all.unadjusted", length(cp.all.unadjusted))))
                  colnames(yy) <- c("y", "grp")
                  yy <- as.data.frame(yy)
                  yy$y <- as.numeric(as.character(yy$y))
                  colnames(yy) <- c("numFeature", "category")
                  m = list(l = 200, r = 5, b = 5, t = 5, pad = 4)
                  p <- plot_ly(yy, x = ~numFeature, color = ~category, type = "box") %>% layout(margin = m)
                  htmlwidgets::saveWidget(p, file.path(output.dir, "boxplot.html"))
                })
                na.pad <- function(x, len) {
                  x[1:len]
                }
                makePaddedDataFrame <- function(l, ...) {
                  maxlen <- max(sapply(l, length))
                  data.frame(lapply(l, na.pad, len = maxlen), ...)
                }
                x <- In.adjusted.not.in.unadjusted
                y <- In.unadjusted.not.in.adjusted
                z <- common
                venn.res <- makePaddedDataFrame(list(x = x, y = y, z = z))
                colnames(venn.res) <- c("adjusted.only", "unadjusted.only", "common")
                write.table(venn.res, file = file.path(output.dir, "results4venn.csv"), quote = FALSE, sep = ",", eol = "\n", na = " ", 
                  dec = ".", row.names = FALSE, col.names = TRUE)
                temp1 <- dplyr::as_data_frame(example.go.adjusted.by.exon[index1, ])
                writeTibble(temp1, file.path(output.dir, "adjustedOnly.csv"))
                temp2 <- dplyr::as_data_frame(example.go.unadjusted[index2, ])
                writeTibble(temp2, file.path(output.dir, "unadjustedOnly.csv"))
            }
            else {
                if (length(In.unadjusted.not.in.adjusted) == 0) {
                  cat("there is no significant gene sets in unadjusted results only\n")
                }
                cat("\n")
                if (length(In.adjusted.not.in.unadjusted) == 0) {
                  cat("there is no significant gene sets in adjusted results only\n")
                }
                cat("\n")
            }
        }
        else {
            cat("The enriched gene sets is less than", n, "\n")
        }
    }
}
gene2cat <- function(gene.name, re) {
    z <- re$genesets
    res <- lapply(z, function(ch) grep(gene.name, ch))
    res2 <- sapply(res, function(x) length(x) > 0)
    gene2cat <- list(re$geneset.names[res2])
    gene2cat
}
gsa.read.gmt <- function(filename) {
    a <- scan(filename, what = list("", ""), sep = "\t", quote = NULL, fill = TRUE, flush = TRUE, multi.line = FALSE)
    geneset.names <- a[1][[1]]
    geneset.descriptions <- a[2][[1]]
    dd <- scan(filename, what = "", sep = "\t", quote = NULL)
    nn <- length(geneset.names)
    n <- length(dd)
    ox <- rep(NA, nn)
    ii <- 1
    for (i in 1:nn) {
        while ((dd[ii] != geneset.names[i]) | (dd[ii + 1] != geneset.descriptions[i])) {
            ii <- ii + 1
        }
        ox[i] <- ii
        ii <- ii + 1
    }
    genesets <- vector("list", nn)
    for (i in 1:(nn - 1)) {
        i1 <- ox[i] + 2
        i2 <- ox[i + 1] - 1
        geneset.descriptions[i] <- dd[ox[i] + 1]
        genesets[[i]] <- dd[i1:i2]
    }
    geneset.descriptions[nn] <- dd[ox[nn] + 1]
    genesets[[nn]] <- dd[(ox[nn] + 2):n]
    out <- list(genesets = genesets, geneset.names = geneset.names, geneset.descriptions = geneset.descriptions)
    class(out) <- "GSA.genesets"
    return(out)
}
gene2cat2 <- function(gmt.input.file) {
    re <- gsa.read.gmt(gmt.input.file)
    gene.name <- unique(do.call(c, re$genesets))
    gene.2.cat <- sapply(gene.name, gene2cat, re)
    names(gene.2.cat) <- gene.name
    gene.2.cat
}
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
reformatdata <- function(re.gene.based) {
    re <- re.gene.based
    no.re.testable.index <- which(re$mostSigDeFeature == "character(0)")
    if (length(no.re.testable.index) > 0) {
        re2 <- re[-no.re.testable.index, ]
    }
    else {
        re2 <- re
    }
    All.gene.id.based.on.sub_feature <- unique(unlist(strsplit(re2$geneID, "\\+")))
    All.gene.id.index <- rep(0, length(All.gene.id.based.on.sub_feature))
    names(All.gene.id.index) <- All.gene.id.based.on.sub_feature
    re3 <- lapply(All.gene.id.based.on.sub_feature, function(u, re2) {
        x <- as.data.frame(re2[grep(u, re2$geneID), ], stringsAsFactors = FALSE)
    }, re2)
    re4 <- do.call(rbind.data.frame, c(re3, stringsAsFactors = FALSE))
    index.geneID <- which(colnames(re4) %in% c("geneID"))
    re5 <- cbind.data.frame(All.gene.id.based.on.sub_feature, re4[, -c(index.geneID)], stringsAsFactors = FALSE)
    colnames(re5)[1] <- "geneID"
    return(re5)
}
heatmap_wPCA = function(Data, g_level = NULL) {
    Data.pca = prcomp(t(Data))
    hmcol <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
    if (is.null(g_level)) {
        type_level = 1:ncol(Data)
        col_level = "black"
        with(data.frame(Data.pca$x), scatter3D(PC1, PC2, PC3, colvar = NULL, type = "h", ticktype = "detailed", bty = "b2", cex = 1, xlab = "PC 1", 
            ylab = "PC 2", zlab = "PC 3", theta = 40, phi = 40, pch = type_level, col = col_level, main = "Principal component analysis"))
        with(data.frame(Data.pca$x), text3D(x = PC1, y = PC2, z = PC3, colnames(Data), col = "black", add = TRUE, colkey = FALSE, cex = 0.5))
    }
    else {
        type_level = 1:ncol(Data)
        TEMP = factor(g_level)
        uniq_label = levels(TEMP)
        levels(TEMP) = hmcol[ceiling(seq(length.out = length(levels(TEMP)), from = 1, to = 256))]
        col_level = as.character(TEMP)
        uniq_col = levels(TEMP)
        Data.pca = prcomp(t(Data))
        with(data.frame(Data.pca$x), scatter3D(PC1, PC2, PC3, colvar = NULL, type = "h", ticktype = "detailed", bty = "b2", cex = 1, xlab = "PC 1", 
            ylab = "PC 2", zlab = "PC 3", theta = 40, phi = 40, pch = type_level, col = col_level, main = "Principal component analysis"))
        legend("topright", legend = uniq_label, pch = type_level, col = uniq_col, cex = 1, inset = c(0.02))
        with(data.frame(Data.pca$x), text3D(x = PC1, y = PC2, z = PC3, colnames(Data), col = "black", add = TRUE, colkey = FALSE, cex = 0.5))
    }
}
writegototable <- function(GO_re, Output_file) {
    dataset2 <- GO_re
    dataset2[sapply(dataset2, is.list)] <- sapply(dataset2[sapply(dataset2, is.list)], function(x) sapply(x, function(y) paste(unlist(y), 
        collapse = ", ")))
    write.table(dataset2, file = Output_file, row.names = FALSE, quote = FALSE, sep = "\t")
}
pathwaysplice <- function(pwf, genome, id, gene2cat, test.cats, go.size.limit, method, repcnt, use.genes.without.cat) {
    if (any(!test.cats %in% c("GO:CC", "GO:BP", "GO:MF", "KEGG"))) {
        stop("Invalid gene_set specified.  Valid categories are GO:CC, GO:BP, GO:MF or KEGG")
    }
    if ((missing(genome) | missing(id))) {
        if (is.null(gene2cat)) {
            stop("You must specify the genome and gene ID format when automatically fetching gene to GO gene_set mappings.")
        }
        genome <- "dummy"
        id <- "dummy"
    }
    if (!any(method %in% c("Wallenius", "Sampling", "Hypergeometric"))) {
        stop("Invalid calculation method selected.  Valid options are Wallenius, Sampling & Hypergeometric.")
    }
    if (!is.null(gene2cat) && (!is.data.frame(gene2cat) & !is.list(gene2cat))) {
        stop("Was expecting a dataframe or a list mapping categories to genes.  Check gene2cat input and try again.")
    }
    pwf <- unfactor(pwf)
    gene2cat <- unfactor(gene2cat)
    if (is.null(gene2cat)) {
        message("Fetching GO annotations...")
        gene2cat <- getGeneSet(rownames(pwf), genome, id, fetch.cats = test.cats)
        cat2gene <- reversemapping(gene2cat)
        gene2cat <- reversemapping(cat2gene)
    }
    else {
        message("Using manually entered categories.")
        if (class(gene2cat) != "list") {
            genecol_sum <- as.numeric(apply(gene2cat, 2, function(u) {
                sum(u %in% rownames(pwf))
            }))
            genecol <- which(genecol_sum != 0)
            if (length(genecol) > 1) {
                genecol <- genecol[order(-genecol_sum)[1]]
                warning(paste("More than one possible gene column found in gene2cat, using the one headed", colnames(gene2cat)[genecol]))
            }
            if (length(genecol) == 0) {
                genecol <- 1
                warning(paste("Gene column could not be identified in gene2cat conclusively, using the one headed", colnames(gene2cat)[genecol]))
            }
            othercol <- 1
            if (genecol == 1) {
                othercol <- 2
            }
            gene2cat <- split(gene2cat[, othercol], gene2cat[, genecol])
            cat2gene <- reversemapping(gene2cat)
            gene2cat <- reversemapping(cat2gene)
        }
        gene2cat <- gene2cat[-which(is.na(names(gene2cat)))]
        if (sum(unique(unlist(gene2cat, use.names = FALSE)) %in% rownames(pwf)) > sum(unique(names(gene2cat)) %in% rownames(pwf))) {
            gene2cat <- reversemapping(gene2cat)
        }
        gene2cat <- gene2cat[names(gene2cat) %in% rownames(pwf)]
        if (length(gene2cat) > 0) {
            cat2gene <- reversemapping(gene2cat)
            gene2cat <- reversemapping(cat2gene)
            cat2gene <- lapply(cat2gene, function(x) {
                unique(x)
            })
            gene2cat <- lapply(gene2cat, function(x) {
                unique(x)
            })
        }
        else {
            cat("There is no match between gene names of gene2pathway input and gene names of the data set under analysis,please change gene2pathway input\n\n")
            return(NA)
        }
    }
    cat2gene <- getGeneSetBySize(cat2gene, go.size.limit)
    if (length(gene2cat) == 0) {
        stop("No gene set is satisfied by the selected size. Change gene set or choose new size.")
    }
    gene2cat <- reversemapping(cat2gene)
    cat2gene <- reversemapping(gene2cat)
    nafrac <- (sum(is.na(pwf$pwf))/nrow(pwf)) * 100
    if (nafrac > 50) {
        warning(paste("Missing length data for ", round(nafrac), "% of genes.  Accuarcy of GO test will be reduced.", sep = ""))
    }
    pwf$pwf[is.na(pwf$pwf)] <- pwf$pwf[match(sort(pwf$bias.data[!is.na(pwf$bias.data)])[ceiling(sum(!is.na(pwf$bias.data))/2)], pwf$bias.data)]
    unknown_go_terms = nrow(pwf) - length(gene2cat)
    if ((!use.genes.without.cat) && unknown_go_terms > 0) {
        message(paste("For", unknown_go_terms, "genes, we could not find any categories. These genes will be excluded."))
        message("To force their use, please run with use_genes_without_cat=TRUE (see documentation).")
        message("This was the default behavior for version 1.15.1 and earlier.")
        pwf = pwf[rownames(pwf) %in% names(gene2cat), ]
    }
    cats <- names(cat2gene)
    DE <- rownames(pwf)[pwf$DEgenes == 1]
    num_de <- length(DE)
    num_genes <- nrow(pwf)
    pvals <- data.frame(gene_set = cats, over_represented_pvalue = NA, under_represented_pvalue = NA, stringsAsFactors = FALSE, numSIGInCat = NA, 
        numInCat = NA)
    if (method == "Sampling") {
        num_DE_mask <- rep(0, length(cats))
        a <- table(unlist(gene2cat[DE], FALSE, FALSE))
        num_DE_mask[match(names(a), cats)] <- as.numeric(a)
        num_DE_mask <- as.integer(num_DE_mask)
        gene2cat <- gene2cat[rownames(pwf)]
        names(gene2cat) <- rownames(pwf)
        message("Running the simulation...")
        lookup <- matrix(0, nrow = repcnt, ncol = length(cats))
        for (i in 1:repcnt) {
            a <- table(as.character(unlist(gene2cat[order(runif(num_genes)^(1/pwf$pwf), decreasing = TRUE)[1:num_de]], FALSE, FALSE)))
            lookup[i, match(names(a), cats)] <- a
            pp(repcnt)
        }
        message("Calculating the p-values...")
        pvals[, 2] <- (colSums(lookup >= outer(rep(1, repcnt), num_DE_mask)) + 1)/(repcnt + 1)
        pvals[, 3] <- (colSums(lookup <= outer(rep(1, repcnt), num_DE_mask)) + 1)/(repcnt + 1)
    }
    if (method == "Wallenius") {
        message("Calculating the p-values...")
        degenesnum <- which(pwf$DEgenes == 1)
        cat2genenum <- relist(match(unlist(cat2gene), rownames(pwf)), cat2gene)
        alpha <- sum(pwf$pwf)
        pvals[, 2:3] <- t(sapply(cat2genenum, function(u) {
            num_de_incat <- sum(degenesnum %in% u)
            num_incat <- length(u)
            avg_weight <- mean(pwf$pwf[u])
            weight <- (avg_weight * (num_genes - num_incat))/(alpha - num_incat * avg_weight)
            if (num_incat == num_genes) {
                weight <- 1
            }
            c(dWNCHypergeo(num_de_incat, num_incat, num_genes - num_incat, num_de, weight) + pWNCHypergeo(num_de_incat, num_incat, num_genes - 
                num_incat, num_de, weight, lower.tail = FALSE), pWNCHypergeo(num_de_incat, num_incat, num_genes - num_incat, num_de, weight))
        }))
    }
    if (method == "Hypergeometric") {
        message("Calculating the p-values...")
        degenesnum <- which(pwf$DEgenes == 1)
        cat2genenum <- relist(match(unlist(cat2gene), rownames(pwf)), cat2gene)
        pvals[, 2:3] <- t(sapply(cat2genenum, function(u) {
            num_de_incat <- sum(degenesnum %in% u)
            num_incat <- length(u)
            c(dhyper(num_de_incat, num_incat, num_genes - num_incat, num_de) + phyper(num_de_incat, num_incat, num_genes - num_incat, num_de, 
                lower.tail = FALSE), phyper(num_de_incat, num_incat, num_genes - num_incat, num_de))
        }))
    }
    degenesnum <- which(pwf$DEgenes == 1)
    cat2genenum <- relist(match(unlist(cat2gene), rownames(pwf)), cat2gene)
    pvals[, 4:5] <- t(sapply(cat2genenum, function(u) {
        c(sum(degenesnum %in% u), length(u))
    }))
    DE_pwf <- rownames(pwf[degenesnum, ])
    pvals.6 <- sapply(cat2gene, function(u, DE_pwf) {
        x <- u[which(u %in% DE_pwf)]
        x
    }, DE_pwf)
    xxx <- match2Genome(genome)
    pvals.6.gene.symbol <- sapply(pvals.6, function(u, xxx) {
        y <- xxx[match(u, as.character(xxx[, 2])), 1]
        y
    }, xxx)
    gene_pwf <- rownames(pwf)
    pvals.7 <- sapply(cat2gene, function(u, gene_pwf) {
        x <- u[which(u %in% gene_pwf)]
        x
    }, gene_pwf)
    pvals.7.gene.symbol <- sapply(pvals.7, function(u, xxx) {
        y <- xxx[match(u, as.character(xxx[, 2])), 1]
        y
    }, xxx)
    pvals.6.df <- list_to_df(pvals.6)
    pvals.6.gene.symbol.df <- list_to_df(pvals.6.gene.symbol)
    pvals.7.df <- list_to_df(pvals.7)
    pvals.7.gene.symbol.df <- list_to_df(pvals.7.gene.symbol)
    dataset2 <- pvals.6.gene.symbol.df
    dataset2[sapply(dataset2, is.list)] <- sapply(dataset2[sapply(dataset2, is.list)], function(x) sapply(x, function(y) paste(unlist(y), 
        collapse = ", ")))
    temp.gene.name <- unique(apply(dataset2[, 2], 1, c))
    temp.gene.name.2 <- unique(gdata::trim(unlist(strsplit(temp.gene.name, split = ","))))
    DE_from_GO <- temp.gene.name.2
    colnames(pvals.6.df) <- c("gene_set", "SIGgene_ensembl")
    colnames(pvals.6.gene.symbol.df) <- c("gene_set", "SIGgene_symbol")
    colnames(pvals.7.df) <- c("gene_set", "All_gene_ensembl")
    colnames(pvals.7.gene.symbol.df) <- c("gene_set", "All_gene_symbol")
    pvals <- pvals[order(pvals$over_represented_pvalue), ]
    if (any(grep("^GO:", pvals$gene_set))) {
        GOnames <- select(GO.db, keys = pvals$gene_set, columns = c("TERM", "ONTOLOGY"))[, 2:3]
        colnames(GOnames) <- tolower(colnames(GOnames))
        colnames(GOnames)[colnames(GOnames) == "term"] <- "description"
        pvals <- cbind(pvals, GOnames)
    }
    re.2 <- merge(pvals, pvals.6.df, by = "gene_set", sort = FALSE)
    re.3 <- merge(re.2, pvals.6.gene.symbol.df, by = "gene_set", sort = FALSE)
    re.4 <- merge(re.3, pvals.7.df, by = "gene_set", sort = FALSE)
    re.5 <- merge(re.4, pvals.7.gene.symbol.df, by = "gene_set", sort = FALSE)
    re.6 <- list(GO = re.5, DE_GO = DE_from_GO, cat2gene = cat2gene)
    return(re.6)
}
getGeneSet <- function(genes, genome, id, fetch.cats = c("GO:CC", "GO:BP", "GO:MF")) {
    if (any(!fetch.cats %in% c("GO:CC", "GO:BP", "GO:MF", "KEGG"))) {
        stop("Invaled gene_set specified.  Categories can only be GO:CC, GO:BP, GO:MF or KEGG")
    }
    orgstring <- as.character(.ORG_PACKAGES[match(gsub("[0-9]+", "", genome), names(.ORG_PACKAGES))])
    if (length(orgstring) != 1) {
        stop("Couldn't grab GO categories automatically.  Please manually specify.")
    }
    library(paste(orgstring, "db", sep = "."), character.only = TRUE)
    coreid <- strsplit(orgstring, "\\.")[[1]][3]
    userid <- as.character(.ID_MAP[match(id, names(.ID_MAP))])
    if (is.na(userid) | (length(userid) != 1)) {
        stop("Couldn't grab GO categories automatically.  Please manually specify.")
    }
    core2cat <- NULL
    if (length(grep("^GO", fetch.cats)) != 0) {
        gomapFunction <- .ORG_GOMAP_FUNCTION[orgstring]
        if (is.na(gomapFunction)) 
            gomapFunction <- .ORG_GOMAP_FUNCTION["default"]
        x <- toTable(get(paste(orgstring, gomapFunction, sep = "")))
        x[!x$Ontology %in% gsub("^GO:", "", fetch.cats), 2] <- "Other"
        core2cat <- x[, 1:2]
        colnames(core2cat) <- c("gene_id", "gene_set")
    }
    if (length(grep("^KEGG", fetch.cats)) != 0) {
        x <- toTable(get(paste(orgstring, "PATH", sep = "")))
        colnames(x) <- c("gene_id", "gene_set")
        if (!is.null(core2cat)) {
            core2cat <- rbind(core2cat, x)
        }
        else {
            core2cat <- x
        }
    }
    if (coreid != userid) {
        user2core <- toTable(get(paste(orgstring, userid, sep = "")))
        user2core <- user2core[user2core[, 1] %in% core2cat[, 1], ]
        list_core2cat <- split(core2cat[, 2], core2cat[, 1])
        list_core2cat <- list_core2cat[match(user2core[, 1], names(list_core2cat))]
        user2cat <- split(unlist(list_core2cat, FALSE, FALSE), rep(user2core[, 2], sapply(list_core2cat, length)))
        user2cat <- sapply(user2cat, unique)
    }
    else {
        user2cat <- split(core2cat[, 2], core2cat[, 1])
        user2cat <- sapply(user2cat, unique)
    }
    user2cat <- lapply(user2cat, function(x) {
        if (length(x) > 1) 
            x = x[x != "Other"]
        x
    })
    names(user2cat) <- toupper(names(user2cat))
    gene2go <- user2cat[toupper(genes)]
    return(gene2go)
}
pp <- function(total, count, i = i) {
    if (missing(count)) {
        count <- evalq(i, envir = parent.frame())
    }
    if (missing(total)) {
        total <- evalq(stop, envir = parent.frame())
    }
    cat(round(100 * (count/total)), "%   \r")
}
outputGoBasedSelection <- function(Re.Go.adjusted.by.exon.SJ) {
    index.select <- which(Re.Go.adjusted.by.exon.SJ[[1]]$numInCat >= 10 & Re.Go.adjusted.by.exon.SJ[[1]]$numInCat <= 300 & Re.Go.adjusted.by.exon.SJ[[1]]$ontology == 
        "BP")
    Re.Go.adjusted.by.exon.SJ.select <- Re.Go.adjusted.by.exon.SJ[[1]][index.select, ]
    Re.Go.adjusted.by.exon.SJ.select <- Re.Go.adjusted.by.exon.SJ.select[, -3]
    temp <- format(Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue, scientific = TRUE, digits = 2)
    Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue <- temp
    rank.value.by.over_represented_pvalue <- rank(as.numeric(Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue), ties.method = "min")
    Re.Go.adjusted.by.exon.SJ.select <- cbind(Re.Go.adjusted.by.exon.SJ.select, rank.value.by.over_represented_pvalue)
    return(Re.Go.adjusted.by.exon.SJ.select)
}
outputCatBasedSelection <- function(Re.Go.adjusted.by.exon.SJ) {
    index.select <- which(Re.Go.adjusted.by.exon.SJ[[1]]$numInCat >= 10 & Re.Go.adjusted.by.exon.SJ[[1]]$numInCat <= 300)
    Re.Go.adjusted.by.exon.SJ.select <- Re.Go.adjusted.by.exon.SJ[[1]][index.select, ]
    Re.Go.adjusted.by.exon.SJ.select <- Re.Go.adjusted.by.exon.SJ.select[, -3]
    temp <- format(Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue, scientific = TRUE, digits = 2)
    Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue <- temp
    rank.value.by.over_represented_pvalue <- rank(as.numeric(Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue), ties.method = "min")
    Re.Go.adjusted.by.exon.SJ.select <- cbind(Re.Go.adjusted.by.exon.SJ.select, rank.value.by.over_represented_pvalue)
    return(Re.Go.adjusted.by.exon.SJ.select)
}
getStaisitcs4Go <- function(GO.wall.DE_interest, mds.11.sample) {
    GO.data <- GO.wall.DE_interest[[1]]
    y <- as.list(GO.data$SIGgene_ensembl)
    re <- lapply(1:length(y), function(u, y, mds.11.sample) {
        yy <- y[[u]]
        y.id <- trim(c(unlist(strsplit(y[[u]], split = ","))))
        if (length(y.id) != 0) {
            yyy <- mean(as.numeric(unlist(mds.11.sample[match(y.id, mds.11.sample$geneID), ]$numFeature)))
        }
        else {
            yyy <- 0
        }
        yyy
    }, y, mds.11.sample)
    re2 <- list_to_df(re)
    GO.data.1 <- cbind(GO.data, re2)
    GO.data.2 <- GO.data.1[, -(dim(GO.data.1)[2] - 1)]
    colnames(GO.data.2)[dim(GO.data.2)[2]] <- "Ave_value_DE"
    cat2gene <- GO.wall.DE_interest[[3]]
    rre <- lapply(1:length(cat2gene), function(u, cat2gene, mds.11.sample) {
        yy <- cat2gene[[u]]
        y.id <- yy
        if (length(y.id) != 0) {
            yyy <- mean(as.numeric(unlist(mds.11.sample[match(y.id, mds.11.sample$geneID), ]$numFeature)), na.rm = TRUE)
        }
        else {
            yyy <- 0
        }
        yyy
    }, cat2gene, mds.11.sample)
    names(rre) <- names(cat2gene)
    rre2 <- list_to_df(rre)
    colnames(rre2) <- c("gene_set", "Ave_value_all_gene")
    GO.data.3 <- merge(GO.data.2, rre2, by = "gene_set", sort = FALSE)
    GO.data.3$Ave_value_DE <- unlist(GO.data.3$Ave_value_DE)
    GO.data.3$Ave_value_all_gene <- unlist(GO.data.3$Ave_value_all_gene)
    re3 <- list(GO.wall.DE_interest = GO.data.3, pwf.DE_interest = GO.wall.DE_interest[[2]])
    return(re3)
}
overlap_ratio <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x, y)))
}
color_scale <- function(c1 = "grey", c2 = "red") {
    pal <- colorRampPalette(c(c1, c2))
    colors <- pal(100)
    return(colors)
}
getIdx <- function(v, MIN, MAX) {
    if (MIN == MAX) {
        return(100)
    }
    intervals <- seq(MIN, MAX, length.out = 100)
    max(which(intervals <= v))
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
        suppressMessages(xxx <- select(edb, keys = entrezid, columns = c("ENSEMBL", "SYMBOL")))
        yyy <- xxx[, c(3, 2, 1)]
    }, mm10 = {
        edb <- org.Mm.eg.db
        entrezid <- keys(edb, keytype = "ENTREZID")
        suppressMessages(xxx <- select(edb, keys = entrezid, columns = c("ENSEMBL", "SYMBOL")))
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
reformatpath <- function(dir.name) {
    CheckOPS <- Sys.info()[["sysname"]]
    if (CheckOPS == "Darwin") {
        temp <- unlist(strsplit(dir.name, split = "\\/"))
        if (!is.na(temp[3] == "H_driver")) {
            if (temp[3] == "H_driver") {
                temp[2] <- "Volumes"
                temp[3] <- "Bioinformatics$"
                dir.name <- do.call("file.path", as.list(temp))
            }
        }
    }
    return(dir.name)
}
getGeneSetBySize <- function(user2cat, go.size.limit) {
    gene2go <- user2cat
    gene2go.select <- lapply(gene2go, function(x) {
        x = x[x != "Other"]
        x
    })
    gene2go.select.1 <- gene2go.select[lapply(gene2go.select, length) > 0]
    lower.size <- go.size.limit[1]
    upper.size <- go.size.limit[2]
    if (is.finite(lower.size) & is.finite(upper.size)) {
        gene2go.select.2 <- gene2go.select.1[lapply(gene2go.select.1, length) > lower.size & lapply(gene2go.select.1, length) <= upper.size]
    }
    else {
        gene2go.select.2 <- gene2go.select.1
    }
    return(gene2go.select.2)
}
makeFeatureTable <- function(jscs, use.multigene.aggregates = FALSE) {
    temp <- fData(jscs)
    temp2 <- temp[which(temp$testable == TRUE), ]
    index.1 <- which(colnames(temp2) %in% c("geneID"))
    index.2 <- which(colnames(temp2) %in% c("countbinID"))
    index.3 <- which(colnames(temp2) %in% c("pvalue"))
    temp3 <- temp2[, c(index.1, index.2, index.3)]
    temp3 <- rapply(temp3, as.character, classes = "factor", how = "replace")
    if (use.multigene.aggregates == FALSE) {
        temp3 <- temp3[-grep("\\+", temp3$geneID), ]
    }
    row.names(temp3) = seq(1, dim(temp3)[1], 1)
    return(temp3)
}
nullpSplice = function(DEgenes, genome, id, bias.data = NULL, plot.fit = TRUE, binsize = "auto") {
    if (!is.null(bias.data) & length(bias.data) != length(DEgenes)) {
        stop("bias.data vector must have the same length as DEgenes vector!")
    }
    bias.data = unfactor(bias.data)
    DEgenes = unfactor(DEgenes)
    if (is.null(bias.data)) {
        bias.data = getlength(names(DEgenes), genome, id)
    }
    pwf = rep(NA, length(DEgenes))
    w = !is.na(bias.data)
    pwf[w] = makespline(bias.data[w], DEgenes[w])
    out = data.frame(DEgenes = DEgenes, bias.data = bias.data, pwf = pwf, stringsAsFactors = FALSE)
    rownames(out) = names(DEgenes)
    if (plot.fit) {
        plotPwfSplice(out, binsize)
    }
    out
}
plotPwfSplice = function(pwf, binsize, pwf_col = 3, pwf_lwd = 2, xlab = "Biased Data in <binsize> gene bins.", ylab = "Proportion of significant genes", 
    ...) {
    w = !is.na(pwf$bias.data)
    o = order(pwf$bias.data[w])
    rang = max(pwf$pwf, na.rm = TRUE) - min(pwf$pwf, na.rm = TRUE)
    if (rang == 0 & binsize == "auto") 
        binsize = 1000
    if (binsize == "auto") {
        binsize = max(1, min(100, floor(sum(w) * 0.08)))
        resid = rang
        oldwarn = options()$warn
        options(warn = -1)
        while (binsize <= floor(sum(w) * 0.1) & resid/rang > 0.001) {
            binsize = binsize + 100
            splitter = ceiling(1:length(pwf$DEgenes[w][o])/binsize)
            de = sapply(split(pwf$DEgenes[w][o], splitter), mean)
            binlen = sapply(split(as.numeric(pwf$bias.data[w][o]), splitter), mean)
            resid = sum((de - approx(pwf$bias.data[w][o], pwf$pwf[w][o], binlen)$y)^2)/length(binlen)
        }
        options(warn = oldwarn)
    }
    else {
        splitter = ceiling(1:length(pwf$DEgenes[w][o])/binsize)
        de = sapply(split(pwf$DEgenes[w][o], splitter), mean)
        binlen = sapply(split(as.numeric(pwf$bias.data[w][o]), splitter), mean)
    }
    xlab = gsub("<binsize>", as.character(binsize), xlab)
    if ("xlab" %in% names(list(...))) {
        if ("ylab" %in% names(list(...))) {
            plot(binlen, de, ...)
        }
        else {
            plot(binlen, de, ylab = ylab, ...)
        }
    }
    else if ("ylab" %in% names(list(...))) {
        plot(binlen, de, xlab = xlab, ...)
    }
    else {
        plot(binlen, de, xlab = xlab, ylab = ylab, ...)
    }
    lines(pwf$bias.data[w][o], pwf$pwf[w][o], col = pwf_col, lwd = pwf_lwd)
}
reformatPathwayOut <- function(pathway.in) {
    res <- dplyr::as_data_frame(pathway.in$GO)
    res
}
writeTibble <- function(tibble.input, output.file.name = tempfile()) {
    if (!dir.exists(dirname(output.file.name))) {
        dir.create(dirname(output.file.name), recursive = TRUE)
    }
    flatten_list = function(x) {
        if (typeof(x) != "list") {
            return(x)
        }
        sapply(x, function(y) paste(y, collapse = " | "))
    }
    tibble.input %>% mutate_all(funs(flatten_list)) %>% write.csv(output.file.name)
}
