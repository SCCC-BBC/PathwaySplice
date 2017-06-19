#' makeGeneTable 
#' 
#' This function obtains genewise p-values,
#' by representing each gene with the smallest p-value among its features, 
#' and then determines genes status as significant or not. 
#'  
#' @param feature.table An \code{featureBasedData} object.
#' @param sig.threshold Significance threshold used to determine whether the gene is significant 
#' or not  
#'
#' @return Returns a genewised table with several variables (columns) 
#' \item{geneID}{Gene ID}
#' \item{geneWisePvalue}{each gene is represented by the smallest p-value among its features}
#' \item{sig.gene}{a gene is significant (1) or not (0)} 
#' \item{mostSigDeFeature}{the most significant gene feature}
#' \item{numFeature}{number of gene features within the gene}
#' 
#' 
#' @export
#'
#' @examples
#' data(featureBasedData)
#' gene.based.table <- makeGeneTable(featureBasedData)
#'   
makeGeneTable <- function(feature.table, sig.threshold = 0.05)
{
    gene.id <- unique(feature.table$geneID)
    
    y <- lapply(gene.id, function(u, feature.table)
    {
        
        x <- as.data.frame(feature.table[which(feature.table$geneID %in% u), 
            ], stringsAsFactors = FALSE)
        
        num.feature <- dim(x)[1]
        
        xx <- x[which.min(x$pvalue), ]
        
        xxx <- cbind.data.frame(xx, num.feature, stringsAsFactors = FALSE)
        
        xxx
        
    }, feature.table)
    
    yy <- do.call(rbind.data.frame, c(y, stringsAsFactors = FALSE))
    
    sig.gene <- ifelse(yy$pvalue < sig.threshold, 1, 0)
    
    z <- cbind.data.frame(yy$geneID, yy$pvalue, sig.gene, yy$countbinID, yy$num.feature, 
        stringsAsFactors = FALSE)
    
    colnames(z) <- c("geneID", "geneWisePvalue", "sig.gene", "mostSigDeFeature", 
        "numFeature")
    
    z <- reformatdata(z)
    
    return(z)
}

#' lrTestBias
#' 
#' This function tests presence of selection bias using logistic regression, and produces a boxplot 
#' that compares distributions of bias factors (e.g. number of exons) for significant genes and non-significant genes. 
#' 
#' @param genewise.table A dataframe with genewise p-value for each gene, returned from \code{makeGeneTable()}
#' @param boxplot.width width of boxplot
#'   
#' @details To determine presentce of selection bias, we fit the logistic regression model 
#' \code{Pr(a gene is significant) ~ number of features within the gene}. 
#' Here features refer to exon bins or splicing junction bins, depending on 
#' how genewise pvalues were obtained in the \code{genewise.table}
#'   
#'   
#' @return Nothing to be returned
#' 
#' @export
#' 
#'  
#'
#' @examples
#' gene.based.table <- makeGeneTable(featureBasedData)
#' lrTestBias(gene.based.table)
#' 
#' 
lrTestBias <- function(genewise.table, boxplot.width = 0.1) 
{
    
    mydata <- genewise.table
    
    n.gene <- dim(mydata)[1]
    
    DE.out <- ifelse(mydata$sig.gene == 1, "Significant genes", "Non-significant genes")
    
    mydata.2 <- cbind(mydata, DE.out)
    
    par(mfrow = c(1, 1))
    
    if (var(as.numeric(unlist(mydata.2$numFeature))) != 0)
    {
        
        mylogit.2 <- glm(DE.out ~ as.numeric(numFeature), data = mydata.2, family = "binomial")
        re <- summary(mylogit.2)
        pvalue <- re$coefficients[2, 4]
        
        p.value <- format.pval(pvalue, eps = .0001, digits = 2)
        
        index.1 <- which(colnames(mydata.2) %in% c("numFeature"))
        index.2 <- which(colnames(mydata.2) %in% c("DE.out"))
        
        temp <- data.frame(mydata.2[, c(index.1, index.2)])
        
        temp$DE.out <- factor(temp$DE.out)
        
        temp$DE.out <- factor(temp$DE.out, levels = levels(temp$DE.out)[c(2, 
            1)])
        
        boxplot(unlist(temp$numFeature) ~ unlist(temp$DE.out), boxwex = boxplot.width, 
            ylab = "Number of features", col = "lightgray", ylim = c(min(temp$numFeature), max(temp$numFeature)))
        
        text(x = 2, y = max(temp$numFeature)-1, labels = c("", paste0("P-value from logistic regression ", 
            p.value)))
    } else
    {
        cat("There are no variations on the number of features\n")
    }
    
}

#' runPathwaySplice
#'
#' This function identifies pathways that are enriched with signficant genes, while accounting for 
#' different number of gene features (e.g. exons) associated with each gene
#' 
#' @param genewise.table data frame returned from function \code{makeGeneTable} 
#' @param genome Genome to be used, options are 'hg19' or 'mm10' 
#' @param id GeneID, options are 'entrezgene' or 'ensembl_gene_id'
#' @param gene2cat Get sets defined by users, can be obtained for example from \code{gmtGene2Cat} function   
#' @param test.cats Default gene sets to be tested if \code{gene2cat} is not defined 
#' @param go.size.limit Size limit of the gene sets to be tested
#' @param method the method used to calculate pathway enrichment p value. 
#'        Options are 'Wallenius', 'Sampling', and 'Hypergeometric' 
#' @param repcnt Number of random samples to be calculated when 'Sampling' is used
#'        ignored unless \code{method='Sampling'}
#' @param use.genes.without.cat Whether genes not mapped to any category tested are included in analysis.
#'        If set to FALSE, genes not mapped to any tested categories are ignored in analysis.
#' @param binsize The number of genes in each gene bin in the bias plot
#' 
#' 
#' @details This function implements the methodology described in Young et al. (2011) to adjust for 
#'          different number of gene features (i.e. counting bins, see Fig 1 in Anders et al. 2012) associated with each gene. In the bias plot, the genes are grouped 
#'          by \code{numFeature} in \code{genewise.table} into gene bins, 
#'          the proportions of signficant genes are then plotted against the gene bins. 
#'              
#' @return runPathwaySplice returns a tibble(data frame) with 11 columns as the following:
#' \item{category}{Name of the gene set (e.g. a pathway or a gene category)} 
#' \item{over_represented_pvalue}{P-vaue for the associated category being over represented among significant genes} 
#' \item{under_represented_pvalue}{P-vaue for the associated category being under represented among significant genes} 
#' \item{numDEInCat}{The number of significant genes in the category} 
#' \item{numInCat}{The total number of genes in the category}                                          
#' \item{term}{The GO term if any of the categories is a GO term} 
#' \item{ontology}{The column for the GO term's ontology if any of the categories is a GO term}
#' \item{DEgene_ID}{The column for the Ensembl gene ID of differentially genes in the category}
#' \item{DEgene_symbol}{The column for the gene symbol of differentially genes in the category}
#' \item{Ave_value_DE}{The column for the average numFeature value of differentially genes in the category}
#' \item{Ave_value_all_gene}{The column for the average numFeature value of total genes in the category}
#'
#' 
#' @references Young MD, Wakefield MJ, Smyth GK, Oshlack A (2011) \emph{Gene ontology analysis for RNA-seq: 
#' accounting for selection bias}. Genome Biology 11:R14
#' 
#' Anders S, Reyes A, Huber W (2012) \emph{Dececting differential usage of exons from RNA-seq data.} 
#' Genome Research 22(10): 2008-2017
#' 
#' @export
#'
#' @examples
#' gene.based.table <- makeGeneTable(featureBasedData)
#' res <- runPathwaySplice(gene.based.table,genome='hg19',id='ensGene',
#'                          test.cats=c('GO:BP'),
#'                          go.size.limit=c(5,30),
#'                          method='Wallenius',binsize=2)
#' 
#'  
runPathwaySplice <- function(genewise.table, genome, id, gene2cat = NULL, test.cats = c("GO:CC", 
    "GO:BP", "GO:MF"), go.size.limit = c(10, 200), method = "Wallenius", repcnt = 2000, 
    use.genes.without.cat = FALSE, binsize = "auto")
    {
    x <- genewise.table$sig.gene
    names(x) <- genewise.table$geneID
    pwf <- nullpSplice(x, genome, id, bias.data = genewise.table$numFeature, plot.fit = TRUE, 
        binsize)
    CatDE <- pathwaysplice(pwf, genome = genome, id = id, gene2cat = gene2cat, 
        test.cats = test.cats, go.size.limit = go.size.limit, method = method, repcnt = repcnt, 
        use.genes.without.cat = use.genes.without.cat)
    res1 <- getStaisitcs4Go(CatDE,genewise.table) 
    res2 <- reformatPathwayOut(res1)
    res2
}

#' enrichmentMap
#'
#' This function draws an enrichment map based on similarities 
#' between gene sets based on Jaccard similarity Coefficient 
#'                                  
#' @param goseqres Object returned from runPathwaySplice
#' @param n The top \emph{n} categories are shown in enrichment map
#' @param fixed If set to FALSE, will invoke tkplot
#' @param vertex.label.font Font size of vertex label
#' @param similarity.threshold Gene sets with Jaccard Coefficient > similarity.threshold will be connected on the enrichment map
#'                
#' @param output.file.dir Output dir for the gene set information file on network
#' @param label.vertex.by.index Which way to be used for labeling vertex on network
#'        
#'        FALSE indicates to label vertex by the name of gene sets
#'        
#'        TRUE indicates to label vertex by the index of gene sets    
#'          
#' @param ... Additional parameter 
#' 
#' @details  The Jaccard similarity coefficient ranges from 0 to 1. JC=0 indicates 
#' there are no overlapping genes between two gene sets, 
#' JC=1 indicates two gene sets are identical. 
#' 
#' @export
#' @return A list for giving edge and vertex information of enrichment map
#' 
#' @author Aimin created this funciton based on enrichMap function in G Yu's DOSE R package
#' 
#' @examples
#' 
#' gene.based.table <- makeGeneTable(featureBasedData)
#' 
#' res <- runPathwaySplice(gene.based.table,genome='hg19',
#'                          id='ensGene',test.cats=c('GO:BP'),
#'                          go.size.limit=c(5,30),
#'                          method='Wallenius')
#' 
#' output.file.dir <- file.path(tempdir(),'OutputEnmapEx')
#' 
#' enmap <- enrichmentMap(res,n=10,similarity.threshold=0,
#'                        output.file.dir = output.file.dir,
#'                        label.vertex.by.index = TRUE)
#'                        
enrichmentMap <- function(goseqres, n = 50, fixed = TRUE, vertex.label.font = 1, 
    similarity.threshold, output.file.dir, label.vertex.by.index = FALSE, ...)
    {
    
    if (!dir.exists(output.file.dir))
    {
        dir.create(output.file.dir,recursive = TRUE)
    }
  
    goseqres <- as.data.frame(goseqres)
  
    GO.name <- goseqres$category
    temp <- goseqres$DEgene_ID
    names(temp) <- GO.name
    x <- goseqres
    geneSets <- temp
    
    y <- as.data.frame(x)
    
    if (any(grep("^GO:", y$category)))
    {
        vertexname <- paste0(y$term, ":", y$numDEInCat)
        
        if (label.vertex.by.index == TRUE)
        {
            vertexname.index <- seq(1, length(vertexname))
            
            output.text <- as.data.frame(cbind(vertexname.index, vertexname))[1:n, 
                ]
            vertexname <- vertexname.index
            colnames(output.text) <- c("index", "name")
            write.table(output.text, file = file.path(output.file.dir, "enrichmap_GO.xls"), 
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
        }
    } else
    {
        vertexname <- paste0(y$category, ":", y$numDEInCat)
        
        if (label.vertex.by.index == TRUE)
        {
            vertexname.index <- seq(1, length(vertexname))
            
            output.text <- as.data.frame(cbind(vertexname.index, vertexname))[1:n, 
                ]
            vertexname <- vertexname.index
            colnames(output.text) <- c("index", "name")
            
            write.table(output.text, file = file.path(output.file.dir, "enrichmap_pathway.xls"), 
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
        }
    }
    
    if (nrow(y) < n)
    {
        n <- nrow(y)
    }
    y <- y[1:n, ]
    
    if (n == 0)
    {
        stop("no enriched term found...")
    } else if (n == 1)
    {
        g <- graph.empty(0, directed = FALSE)
        g <- add_vertices(g, nv = 1)
        
        igraph::V(g)$name <- vertexname
        igraph::V(g)$color <- "red"
        
    } else
    {
        pvalue <- as.numeric(y$over_represented_pvalue)
        
        id <- y[, 1]
        geneSets <- geneSets[id]
        
        n <- nrow(y)  #
        w <- matrix(NA, nrow = n, ncol = n)
        colnames(w) <- rownames(w) <- vertexname[1:n]
        
        for (i in 1:n)
        {
            for (j in i:n)
            {
                w[i, j] <- overlap_ratio(geneSets[id[i]], geneSets[id[j]])
            }
        }
        
        wd <- melt(w)
        wd <- wd[wd[, 1] != wd[, 2], ]
        wd <- wd[!is.na(wd[, 3]), ]
        g <- graph.data.frame(wd[, -3], directed = FALSE)
        igraph::E(g)$width <- sqrt(wd[, 3] * 5)
        
        g <- delete.edges(g, igraph::E(g)[wd[, 3] < similarity.threshold])
        
        idx <- unlist(sapply(igraph::V(g)$name, function(x) which(x == vertexname[1:n])))
        
        cols <- color_scale("red", "#E5C494")
        
        igraph::V(g)$color <- cols[sapply(pvalue, getIdx, min(pvalue), max(pvalue))]
        
        Edata <- as.data.frame(get.edgelist(g))
        Edata$edgewidth <- igraph::E(g)$width
        Vdata <- data.frame(pathway = igraph::V(g)$name, color = igraph::V(g)$color)
        map_data <- list(edge_data = Edata, vertex_data = Vdata)
        
        cnt <- as.integer(y$numDEInCat)
        
        names(cnt) <- vertexname[1:n]
        
        cnt2 <- cnt[igraph::V(g)$name]
        
        igraph::V(g)$size <- cnt2/sum(cnt2) * 100
    }
    
    netplot(g, vertex.label.font = vertex.label.font, vertex.label.color = "black", 
        fixed = fixed, ...)
    
    invisible(g)
    
    re2 <- map_data
    return(re2)
}

#' gmtGene2Cat
#'
#' Obtains all pathways associated with a set of genes 
#' 
#' @param dir.name Directory where the gene sets file (in GMT format) is located 
#' @param pathway.file Input file for the gene sets in GMT format
#' @param location.type Location of the gene set file. Valid options are "local" or "url"
#' @param gene.anno.file Gene annotation file supplied as a file 
#' @param genomeID Genome ('mm10','hg19' or 'hg38') to be used
#'
#' @details This function reads a gene set file in GMT format (http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29), and returns a list with its name
#' being a gene id, and each element of the list being the pathways associated with the gene
#'
#' @return A list where each entry is named by a gene and contains a vector of all
#'         the pathways associated with the gene
#'
#' @export
#'
#' @examples
#' 
#' dir.name <- system.file('extdata', package='PathwaySplice')
#' canonical.pathway.file <- '10.cp.gmt.txt'
#' cpp <- gmtGene2Cat(dir.name,canonical.pathway.file,'local',genomeID='hg19')
#' 
gmtGene2Cat <- function(dir.name, pathway.file, location.type, gene.anno.file = NULL, 
    genomeID = c("mm10", "hg19", "hg38"))
    {
    
    gmt_input_file <- file.path(dir.name, pathway.file)
    
    gene.2.cat.gmt <- gene2cat2(gmt_input_file, location.type)
    
    names.gene.gmt <- as.data.frame(names(gene.2.cat.gmt))
    colnames(names.gene.gmt) <- "gene_id"
    
    if (!is.null(gene.anno.file))
    {
        gene.anno.dir <- dirname(gene_anno_file)
        gene.annno.dir <- reformatpath(gene.anno.dir)
        file.name <- basename(gene_anno_file)
        
        gene_anno_file <- file.path(dir.name, file.name)
        
        gene.id.conversion <- read.csv(gene_anno_file)
    } else
    {
        gene.id.conversion <- match.arg(genomeID)
    }
    
    xxx <- match2Genome(gene.id.conversion)
    
    names.gene.gmt.2 <- match(names.gene.gmt$gene_id, xxx[, 1])
    
    gene.id.conversion.2 <- xxx[names.gene.gmt.2, ]
    
    gene.2.cat.gmt.2 <- gene.2.cat.gmt
    names(gene.2.cat.gmt.2) <- gene.id.conversion.2[, 2]
    gene.2.cat.gmt.2
    
}

#' compareResults
#'
#' @param n.go Number of gene sets
#' @param adjusted Adjusted result 
#' @param unadjusted Unadjusted result 
#' @param output.dir Path for output( including venn,boxplot and two tables)
#' @param type.boxplot Get boxplot for 5 categories or 3 categories
#' 
#'        6 categories: 'All.adjusted','All.unadjusted',
#'                      'Top25.adjusted','Top25.unadjusted',
#'                      'In_ad_not_un','In_un_not_ad'
#'                      
#'        3 categories: 'All','Top25.adjusted','Top25.unadjusted
#'        
#' @param In.ad.not.un.file File name for outputing adjused but not
#'        in unadjusted when using the selected gene sets     
#' @param In.un.not.ad.file File name for outputing unadjused but not
#'        in adjusted when using the selected gene sets  
#'
#' @return Nothing to be returned
#' 
#' @examples
#' 
#' dir.name <- system.file('extdata', package='PathwaySplice')
#' 
#' canonical.pathway.file <- '10.cp.gmt.txt'
#' 
#' cpp <- gmtGene2Cat(dir.name,canonical.pathway.file,
#'                    'local',genomeID='hg19')
#'                    
#' gene.based.table <- makeGeneTable(featureBasedData)
#' 
#' res1 <- runPathwaySplice(gene.based.table,genome='hg19',
#'                          id='ensGene',gene2cat=cpp,go.size.limit = c(2, 200),
#'                          method='Wallenius')
#' 
#' res2 <- runPathwaySplice(gene.based.table,genome='hg19',
#'                          id='ensGene',gene2cat=cpp,go.size.limit = c(2, 200),
#'                          method='Hypergeometric')
#' 
#' output.dir <- file.path(tempdir(),'OutputPostAnalysis')
#' 
#' output.file.name.1 <- 'In_ad_not_un.xls'
#' output.file.name.2 <- 'In_un_not_ad.xls'
#' 
#' compareResults(4,res1,res2,output.dir,
#'                       type.boxplot='Only3',
#'                       output.file.name.1,output.file.name.2)
#' @export
compareResults <- function(n.go, adjusted,unadjusted,output.dir, 
    type.boxplot = c("All", "Only3"), In.ad.not.un.file, In.un.not.ad.file)
    {
    
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir,recursive = TRUE)
    }
    
    n <- n.go
    
    example.go.adjusted.by.exon <- as.data.frame(adjusted)
    example.go.unadjusted <- as.data.frame(unadjusted)
    
    if (is.na(example.go.adjusted.by.exon) || is.na(example.go.unadjusted))
    {
        cat("One of results is empty\n\n")
        return()
    } else
    {

        if (dim(example.go.adjusted.by.exon)[1] >= n && dim(example.go.unadjusted)[1] >= 
            n)
            {
            
            adjusted <- example.go.adjusted.by.exon[1:n, 1]
            unadjusted <- example.go.unadjusted[1:n, 1]
            
            re <- list(adjusted = adjusted, unadjusted = unadjusted)
            
            venn.plot <- venn.diagram(x = re[c(1, 2)], filename = file.path(output.dir, 
                paste0(names(re)[1], "_", names(re)[2], "_overlap_venn.tiff")), 
                height = 3000, width = 3500, resolution = 1000, col = "black", 
                lty = "dotted", lwd = 1, fill = c("red", "blue"), alpha = 0.5, 
                label.col = c(rep("black", 3)), cex = 0.5, fontfamily = "serif", 
                fontface = "bold", cat.col = c("red", "blue"), cat.cex = 0.5, 
                cat.pos = 0.5, cat.dist = 0.05, cat.fontfamily = "serif")
            
            # boxplot
            common <- intersect(unadjusted, adjusted)
            
            In.unadjusted.not.in.adjusted <- setdiff(unadjusted, common)
            In.adjusted.not.in.unadjusted <- setdiff(adjusted, common)
            
            if (length(In.unadjusted.not.in.adjusted) != 0 && length(In.adjusted.not.in.unadjusted) != 
                0)
                {
                index1 <- match(In.adjusted.not.in.unadjusted, example.go.adjusted.by.exon$category)
                In.ad.not.un <- example.go.adjusted.by.exon[index1, ]$Ave_value_all_gene
                
                yy <- cbind(example.go.unadjusted[index1, ]$rank.value.by.over_represented_pvalue, 
                  example.go.adjusted.by.exon[index1, ]$rank.value.by.over_represented_pvalue)
                
                
                index2 <- match(In.unadjusted.not.in.adjusted, example.go.unadjusted$category)
                In.un.not.ad <- example.go.unadjusted[index2, ]$Ave_value_all_gene
                
                yyy <- cbind(example.go.unadjusted[index2, ]$rank.value.by.over_represented_pvalue, 
                  example.go.adjusted.by.exon[index2, ]$rank.value.by.over_represented_pvalue)
                
                xx <- cbind(unlist(In.ad.not.un), unlist(In.un.not.ad))
                
                colnames(xx) <- c("In.ad.not.un", "In.un.not.ad")
                
                cp.topN.adjusted <- unlist(example.go.adjusted.by.exon[1:n, 
                  ]$Ave_value_all_gene)
                cp.topN.unadjusted <- unlist(example.go.unadjusted[1:n, 
                  ]$Ave_value_all_gene)
                
                cp.all.adjusted <- unlist(example.go.adjusted.by.exon$Ave_value_all_gene)
                cp.all.unadjusted <- unlist(example.go.unadjusted$Ave_value_all_gene)
                
                type.boxplot <- match.arg(type.boxplot)
                
                switch(type.boxplot, Only3 = {
                  yy <- rbind(cbind(cp.topN.adjusted, rep(paste0("Adjusted_",n.go), length(cp.topN.adjusted))), 
                    cbind(cp.topN.unadjusted, rep(paste0("Unadjusted_",n.go), length(cp.topN.unadjusted))), 
                    cbind(cp.all.unadjusted, rep("All", length(cp.all.unadjusted))))
                  colnames(yy) <- c("y", "grp")
                  yy <- as.data.frame(yy)
                  yy$grp <- factor(yy$grp)
                  yy$grp <- factor(yy$grp, levels = levels(yy$grp)[c(2, 1, 3)])
                  png(file.path(output.dir, "boxplot.png"))
                  boxplot(as.numeric(as.character(y)) ~ grp, data = yy)
                  dev.off()
                }, {
                  yy <- rbind(cbind(xx[, 1], rep("In.ad.not.un", length(xx[, 
                    1]))), cbind(xx[, 2], rep("In.un.not.ad", length(xx[, 2]))), 
                    cbind(cp.topN.adjusted, rep(paste0("cp.top.adjusted.",n.go), length(cp.topN.adjusted))), 
                    cbind(cp.topN.unadjusted, rep(paste0("cp.top.unadjusted.",n.go), 
                      length(cp.topN.unadjusted))), cbind(cp.all.adjusted, 
                      rep("cp.all.adjusted", length(cp.all.adjusted))), cbind(cp.all.unadjusted, 
                      rep("cp.all.unadjusted", length(cp.all.unadjusted))))
                  colnames(yy) <- c("y", "grp")
                  yy <- as.data.frame(yy)
                  png(file.path(output.dir, "boxplot.png"))
                  boxplot(as.numeric(as.character(y)) ~ grp, data = yy)
                  dev.off()
                })
                
                Output_file <- file.path(output.dir, In.ad.not.un.file)
                writegototable(example.go.adjusted.by.exon[index1, ], Output_file)
                
                Output_file <- file.path(output.dir, In.un.not.ad.file)
                writegototable(example.go.unadjusted[index2, ], Output_file)
                
            } else
            {
                
                if (length(In.unadjusted.not.in.adjusted) == 0)
                {
                  cat("there is no gene sets in unadjusted resutls but not in adjusted resutls\n")
                }
                
                cat("\n")
                
                if (length(In.adjusted.not.in.unadjusted) == 0)
                {
                  cat("there is no gene sets in adjusted resutls but not in unadjusted resutls\n")
                }
                
                cat("\n")
                
            }
        } else
        {
            
            cat("The enriched gene sets is less than", n, "\n")
            
        }
    }
}

#' Title
#'
#' @param libname
#' @param pkgname
#'
#' @return Nothing to be returned
#' @export
#'
#' @examples

.onAttach <- function(libname, pkgname)
{
    if (.Platform$OS.type == "windows" && .Platform$GUI == "Rgui")
    {
        winMenuAddItem("Vignettes", "PathwaySplice", "shell.exec(system.file(\"doc\",\"PathwaySplice.pdf\",package=\"PathwaySplice\"))")
    }
}

# These two variables are required for automatic fetching of categories to
# function.  Their purpose is to take the UCSC genome and gene ID values
# given when looking up length data and convert them to the names used for
# the same organism and gene identifier in the organism packages.

# Mappings that are primarily required by getgo, the purpose of this is to
# convert the UCSC genome IDs, to the bioconductor organism names, e.g.
# 'mm'->'org.Mm.'
.ORG_PACKAGES = paste("org.", c("Ag.eg", "At.tair", "Bt.eg", "Ce.eg", "Cf.eg", 
    "Dm.eg", "Dr.eg", "EcK12.eg", "EcSakai.eg", "Gg.eg", "Hs.eg", "Mm.eg", "Mmu.eg", 
    "Pf.plasmo", "Pt.eg", "Rn.eg", "Sc.sgd", "Ss.eg", "Xl.eg"), sep = "")
names(.ORG_PACKAGES) = c("anoGam", "Arabidopsis", "bosTau", "ce", "canFam", 
    "dm", "danRer", "E. coli K12", "E. coli Sakai", "galGal", "hg", "mm", "rheMac", 
    "Malaria", "panTro", "rn", "sacCer", "susScr", "xenTro")

# These are the only formats supported by getgo at the moment, the purpose
# is to convert the USCC gene ID formats, to the shorthand used by the
# bioconductor organism packages, .e.g. 'refGene'->'ENSEMBL'
.ID_MAP = c("eg", "eg", "ENSEMBL", "SYMBOL", "sgd", "plasmo", "tair")
names(.ID_MAP) = c("knownGene", "refGene", "ensGene", "geneSymbol", "sgd", "plasmo", 
    "tair")

# Below are the exceptions to the function name for gene to go term mappings
.ORG_GOMAP_FUNCTION = c("GO2ALLEGS", "GO2ALLTAIRS", "GO2ALLORFS", "GO2ALLORFS")
names(.ORG_GOMAP_FUNCTION) = c("default", "org.At.tair", "org.Pf.plasmo", "org.Sc.sgd")

# TxDb Length databases
.TXDB_ORGS = c("ce6", "dm3", "hg18", "hg19", "hg38", "mm10", "mm9", "rn4", "rn5", 
    "sacCer2", "sacCer3")

# Utility functions for PathwaySplice
gene2cat <- function(gene.name, re)
{
    z <- re$genesets
    res <- lapply(z, function(ch) grep(gene.name, ch))
    res2 <- sapply(res, function(x) length(x) > 0)
    gene2cat <- list(re$geneset.names[res2])
    gene2cat
}

gsa.read.gmt <- function(filename, file.location.option=c("url","local"))
{
  
  file.location.option <- match.arg(file.location.option)
  
  switch (file.location.option,
          
          url = {
            filename <- filename
          },
          {
            dir.name <- dirname(filename)
            dir.name <- reformatpath(dir.name)
            file.name <- basename(filename)
            filename <- file.path(dir.name, file.name)
          }
  )
  
    a <- scan(filename, what = list("", ""), sep = "\t", quote = NULL, fill = TRUE, 
        flush = TRUE, multi.line = FALSE)
    geneset.names <- a[1][[1]]
    geneset.descriptions <- a[2][[1]]
    dd <- scan(filename, what = "", sep = "\t", quote = NULL)
    nn <- length(geneset.names)
    n <- length(dd)
    ox <- rep(NA, nn)
    ii <- 1
    for (i in 1:nn)
    {
        while ((dd[ii] != geneset.names[i]) | (dd[ii + 1] != geneset.descriptions[i]))
        {
            ii <- ii + 1
        }
        ox[i] <- ii
        ii <- ii + 1
    }
    genesets <- vector("list", nn)
    for (i in 1:(nn - 1))
    {
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

gene2cat2 <- function(gmt.input.file, location.type)
{
    
    re <- gsa.read.gmt(gmt.input.file, location.type)
    gene.name <- unique(do.call(c, re$genesets))
    gene.2.cat <- sapply(gene.name, gene2cat, re)
    names(gene.2.cat) <- gene.name
    gene.2.cat
    
}

list_to_df <- function(list_for_df)
{
    list_for_df <- as.list(list_for_df)
    
    nm <- names(list_for_df)
    if (is.null(nm)) 
        nm <- seq_along(list_for_df)
    
    df <- data.frame(name = nm, stringsAsFactors = FALSE)
    df$value <- unname(list_for_df)
    df
}

reversemapping <- function(map)
{
    tmp <- unlist(map, use.names = FALSE)
    names(tmp) <- rep(names(map), times = as.numeric(summary(map)[, 1]))
    return(split(names(tmp), as.vector(tmp)))
}

reformatdata <- function(re.gene.based)
{
    # re <- pData(re.PJ.gene.based)
    re <- re.gene.based
    no.re.testable.index <- which(re$mostSigDeFeature == "character(0)")
    
    if (length(no.re.testable.index) > 0)
    {
        re2 <- re[-no.re.testable.index, ]
    } else
    {
        re2 <- re
    }
    
    All.gene.id.based.on.sub_feature <- unique(unlist(strsplit(re2$geneID, "\\+")))
    
    All.gene.id.index <- rep(0, length(All.gene.id.based.on.sub_feature))
    names(All.gene.id.index) <- All.gene.id.based.on.sub_feature
    
    re3 <- lapply(All.gene.id.based.on.sub_feature, function(u, re2)
    {
        x <- as.data.frame(re2[grep(u, re2$geneID), ], stringsAsFactors = FALSE)
    }, re2)
    
    re4 <- do.call(rbind.data.frame, c(re3, stringsAsFactors = FALSE))
    
    index.geneID <- which(colnames(re4) %in% c("geneID"))
    re5 <- cbind.data.frame(All.gene.id.based.on.sub_feature, re4[, -c(index.geneID)], 
        stringsAsFactors = FALSE)
    colnames(re5)[1] <- "geneID"
    
    return(re5)
    
}

heatmap_wPCA = function(Data, g_level = NULL)
{
    
    Data.pca = prcomp(t(Data))
    hmcol <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
    
    if (is.null(g_level))
    {
        type_level = 1:ncol(Data)
        col_level = "black"
        
        with(data.frame(Data.pca$x), scatter3D(PC1, PC2, PC3, colvar = NULL, 
            type = "h", ticktype = "detailed", bty = "b2", cex = 1, xlab = "PC 1", 
            ylab = "PC 2", zlab = "PC 3", theta = 40, phi = 40, pch = type_level, 
            col = col_level, main = "Principal component analysis"))
        
        
        
        with(data.frame(Data.pca$x), text3D(x = PC1, y = PC2, z = PC3, colnames(Data), 
            col = "black", add = TRUE, colkey = FALSE, cex = 0.5))
    } else
    {
        type_level = 1:ncol(Data)
        TEMP = factor(g_level)
        uniq_label = levels(TEMP)
        levels(TEMP) = hmcol[ceiling(seq(length.out = length(levels(TEMP)), 
            from = 1, to = 256))]
        col_level = as.character(TEMP)
        uniq_col = levels(TEMP)
        
        Data.pca = prcomp(t(Data))
        with(data.frame(Data.pca$x), scatter3D(PC1, PC2, PC3, colvar = NULL, 
            type = "h", ticktype = "detailed", bty = "b2", cex = 1, xlab = "PC 1", 
            ylab = "PC 2", zlab = "PC 3", theta = 40, phi = 40, pch = type_level, 
            col = col_level, main = "Principal component analysis"))
        
        legend("topright", legend = uniq_label, pch = type_level, col = uniq_col, 
            cex = 1, inset = c(0.02))
        
        with(data.frame(Data.pca$x), text3D(x = PC1, y = PC2, z = PC3, colnames(Data), 
            col = "black", add = TRUE, colkey = FALSE, cex = 0.5))
    }
}

writegototable <- function(GO_re, Output_file)
{
    dataset2 <- GO_re
    dataset2[sapply(dataset2, is.list)] <- sapply(dataset2[sapply(dataset2, 
        is.list)], function(x) sapply(x, function(y) paste(unlist(y), collapse = ", ")))
    
    write.table(dataset2, file = Output_file, row.names = FALSE, quote = FALSE, 
        sep = "\t")
}

pathwaysplice <- function(pwf, genome, id, gene2cat, test.cats, go.size.limit, 
    method, repcnt, use.genes.without.cat)
    {
    ################# Input pre-processing and validation ################### Do some validation
    ################# of input variables
    if (any(!test.cats %in% c("GO:CC", "GO:BP", "GO:MF", "KEGG")))
    {
        stop("Invalid category specified.  Valid categories are GO:CC, GO:BP, GO:MF or KEGG")
    }
    if ((missing(genome) | missing(id)))
    {
        if (is.null(gene2cat))
        {
            stop("You must specify the genome and gene ID format when automatically fetching gene to GO category mappings.")
        }
        # If we're using user specified mappings, this obviously isn't a problem
        genome <- "dummy"
        id <- "dummy"
    }
    if (!any(method %in% c("Wallenius", "Sampling", "Hypergeometric")))
    {
        stop("Invalid calculation method selected.  Valid options are Wallenius, Sampling & Hypergeometric.")
    }
    if (!is.null(gene2cat) && (!is.data.frame(gene2cat) & !is.list(gene2cat)))
    {
        stop("Was expecting a dataframe or a list mapping categories to genes.  Check gene2cat input and try again.")
    }
    
    # Factors are evil
    pwf <- unfactor(pwf)
    gene2cat <- unfactor(gene2cat)
    
    ###################### Data fetching and processing ########################
    if (is.null(gene2cat))
    {
        # When we fetch the data using getgo it will be in the list format
        message("Fetching GO annotations...")
        gene2cat <- getGeneSet(rownames(pwf), genome, id, fetch.cats = test.cats)
        # names(gene2cat) <- rownames(pwf)
        
        # cat('OK') Do the two rebuilds to remove any nulls
        cat2gene <- reversemapping(gene2cat)
        gene2cat <- reversemapping(cat2gene)
        
        # print(cat2gene) print(gene2cat)
        
    } else
    {
        # The gene2cat input accepts a number of formats, we need to check each of
        # them in term
        message("Using manually entered categories.")
        # The options are a flat mapping (that is a data frame or matrix) or a list,
        # where the list can be either gene->categories or category->genes
        if (class(gene2cat) != "list")
        {
            # it's not a list so it must be a data.frame, work out which column contains
            # the genes
            genecol_sum <- as.numeric(apply(gene2cat, 2, function(u)
            {
                sum(u %in% rownames(pwf))
            }))
            genecol <- which(genecol_sum != 0)
            if (length(genecol) > 1)
            {
                genecol <- genecol[order(-genecol_sum)[1]]
                warning(paste("More than one possible gene column found in gene2cat, using the one headed", 
                  colnames(gene2cat)[genecol]))
            }
            if (length(genecol) == 0)
            {
                genecol <- 1
                warning(paste("Gene column could not be identified in gene2cat conclusively, using the one headed", 
                  colnames(gene2cat)[genecol]))
            }
            othercol <- 1
            if (genecol == 1)
            {
                othercol <- 2
            }
            # Now put it into our delicious listy format
            gene2cat <- split(gene2cat[, othercol], gene2cat[, genecol])
            # Do the appropriate builds
            cat2gene <- reversemapping(gene2cat)
            gene2cat <- reversemapping(cat2gene)
            
        }
        
        gene2cat <- gene2cat[-which(is.na(names(gene2cat)))]
        
        # !!!! The following conditional has been flagged as a potential issue when
        # using certain types of input where the category names are the same as gene
        # names (which seems like something you should avoid anyway...).  Leave it
        # for now !!!!  We're now garunteed to have a list (unless the user screwed
        # up the input) but it could be category->genes rather than the
        # gene->categories that we want.
        if (sum(unique(unlist(gene2cat, use.names = FALSE)) %in% rownames(pwf)) > 
            sum(unique(names(gene2cat)) %in% rownames(pwf)))
            {
            gene2cat <- reversemapping(gene2cat)
        }
        # Alright, we're garunteed a list going in the direction we want now.  Throw
        # out genes which we will not use
        gene2cat <- gene2cat[names(gene2cat) %in% rownames(pwf)]
        
        if (length(gene2cat) > 0)
        {
            
            # Rebuild because it's a fun thing to do
            cat2gene <- reversemapping(gene2cat)
            gene2cat <- reversemapping(cat2gene)
            
            ## make sure we remove duplicate entries .. e.g. see
            ## http://permalink.gmane.org/gmane.science.biology.informatics.conductor/46876
            cat2gene <- lapply(cat2gene, function(x)
            {
                unique(x)
            })
            gene2cat <- lapply(gene2cat, function(x)
            {
                unique(x)
            })
        } else
        {
            
            cat("There is no match between gene names of gene2pathway input and gene names of the data set under analysis,please change gene2pathway input\n\n")
            
            return(NA)
            
        }
        
    }
    
    # Add option to choose gene set by its size
    gene2cat <- getGeneSetBySize(gene2cat, go.size.limit)
    
    if(length(gene2cat)==0){
      stop("No gene set is satisfied by the selected size. Change gene set or choose new size.")
    }
    
    cat2gene <- reversemapping(gene2cat)
    gene2cat <- reversemapping(cat2gene)
    
    nafrac <- (sum(is.na(pwf$pwf))/nrow(pwf)) * 100
    if (nafrac > 50)
    {
        warning(paste("Missing length data for ", round(nafrac), "% of genes.  Accuarcy of GO test will be reduced.", 
            sep = ""))
    }
    # Give the genes with unknown length the weight used by the median gene (not
    # the median weighting!)
    pwf$pwf[is.na(pwf$pwf)] <- pwf$pwf[match(sort(pwf$bias.data[!is.na(pwf$bias.data)])[ceiling(sum(!is.na(pwf$bias.data))/2)], 
        pwf$bias.data)]
    
    
    
    # Remove all the genes with unknown GOterms
    unknown_go_terms = nrow(pwf) - length(gene2cat)
    if ((!use.genes.without.cat) && unknown_go_terms > 0)
    {
        message(paste("For", unknown_go_terms, "genes, we could not find any categories. These genes will be excluded."))
        message("To force their use, please run with use_genes_without_cat=TRUE (see documentation).")
        message("This was the default behavior for version 1.15.1 and earlier.")
        pwf = pwf[rownames(pwf) %in% names(gene2cat), ]
    }
    
    # A few variables are always useful so calculate them
    cats <- names(cat2gene)
    DE <- rownames(pwf)[pwf$DEgenes == 1]
    num_de <- length(DE)
    num_genes <- nrow(pwf)
    pvals <- data.frame(category = cats, over_represented_pvalue = NA, under_represented_pvalue = NA, 
        stringsAsFactors = FALSE, numDEInCat = NA, numInCat = NA)
    if (method == "Sampling")
    {
        # We need to know the number of DE genes in each category, make this as a
        # mask that we can use later...
        num_DE_mask <- rep(0, length(cats))
        a <- table(unlist(gene2cat[DE], FALSE, FALSE))
        
        num_DE_mask[match(names(a), cats)] <- as.numeric(a)
        num_DE_mask <- as.integer(num_DE_mask)
        # We have to ensure that genes not associated with a category are included
        # in the simulation, to do this they need an empty entry in the gene2cat
        # list
        gene2cat <- gene2cat[rownames(pwf)]
        names(gene2cat) <- rownames(pwf)
        message("Running the simulation...")
        # Now do the actual simulating
        lookup <- matrix(0, nrow = repcnt, ncol = length(cats))
        for (i in 1:repcnt)
        {
            # A more efficient way of doing weighted random sampling without replacment
            # than the built in function The order(runif...)[1:n] bit picks n genes at
            # random, weighting them by the PWF The table(as.character(unlist(...))) bit
            # then counts the number of times this random set occured in each category
            a <- table(as.character(unlist(gene2cat[order(runif(num_genes)^(1/pwf$pwf), 
                decreasing = TRUE)[1:num_de]], FALSE, FALSE)))
            lookup[i, match(names(a), cats)] <- a
            pp(repcnt)
        }
        message("Calculating the p-values...")
        # The only advantage of the loop is it uses less memory...  for(i in
        # 1:length(cats)){
        # pvals[i,2:3]=c((sum(lookup[,i]>=num_DE_mask[i])+1)/(repcnt+1),(sum(lookup[,i]<=num_DE_mask[i])+1)/(repcnt+1))
        # pp(length(cats)) }
        pvals[, 2] <- (colSums(lookup >= outer(rep(1, repcnt), num_DE_mask)) + 
            1)/(repcnt + 1)
        pvals[, 3] <- (colSums(lookup <= outer(rep(1, repcnt), num_DE_mask)) + 
            1)/(repcnt + 1)
    }
    if (method == "Wallenius")
    {
        message("Calculating the p-values...")
        # All these things are just to make stuff run faster, mostly because
        # comparison of integers is faster than string comparison
        degenesnum <- which(pwf$DEgenes == 1)
        # Turn all genes into a reference to the pwf object
        
        cat2genenum <- relist(match(unlist(cat2gene), rownames(pwf)), cat2gene)
        # This value is used in every calculation, by storing it we need only
        # calculate it once
        alpha <- sum(pwf$pwf)
        
        # Each category will have a different weighting so needs its own test
        pvals[, 2:3] <- t(sapply(cat2genenum, function(u)
        {
            # The number of DE genes in this category
            num_de_incat <- sum(degenesnum %in% u)
            
            # The total number of genes in this category
            num_incat <- length(u)
            
            # This is just a quick way of calculating weight=avg(PWF within
            # category)/avg(PWF outside of category)
            avg_weight <- mean(pwf$pwf[u])
            weight <- (avg_weight * (num_genes - num_incat))/(alpha - num_incat * 
                avg_weight)
            if (num_incat == num_genes) 
                {
                  weight <- 1
                }  #case for the root GO terms
            
            # Now calculate the sum of the tails of the Wallenius distribution (the
            # p-values)
            
            c(dWNCHypergeo(num_de_incat, num_incat, num_genes - num_incat, num_de, 
                weight) + pWNCHypergeo(num_de_incat, num_incat, num_genes - 
                num_incat, num_de, weight, lower.tail = FALSE), pWNCHypergeo(num_de_incat, 
                num_incat, num_genes - num_incat, num_de, weight))
        }))
    }
    if (method == "Hypergeometric")
    {
        message("Calculating the p-values...")
        # All these things are just to make stuff run faster, mostly because
        # comparison of integers is faster than string comparison
        degenesnum <- which(pwf$DEgenes == 1)
        # Turn all genes into a reference to the pwf object
        cat2genenum <- relist(match(unlist(cat2gene), rownames(pwf)), cat2gene)
        # Simple hypergeometric test, one category at a time
        pvals[, 2:3] <- t(sapply(cat2genenum, function(u)
        {
            # The number of DE genes in this category
            num_de_incat <- sum(degenesnum %in% u)
            # The total number of genes in this category
            num_incat <- length(u)
            # Calculate the sum of the tails of the hypergeometric distribution (the
            # p-values)
            c(dhyper(num_de_incat, num_incat, num_genes - num_incat, num_de) + 
                phyper(num_de_incat, num_incat, num_genes - num_incat, num_de, 
                  lower.tail = FALSE), phyper(num_de_incat, num_incat, num_genes - 
                num_incat, num_de))
        }))
    }
    
    # Populate the count columns...
    degenesnum <- which(pwf$DEgenes == 1)
    cat2genenum <- relist(match(unlist(cat2gene), rownames(pwf)), cat2gene)
    pvals[, 4:5] <- t(sapply(cat2genenum, function(u)
    {
        c(sum(degenesnum %in% u), length(u))
    }))
    
    DE_pwf <- rownames(pwf[degenesnum, ])
    
    pvals.6 <- sapply(cat2gene, function(u, DE_pwf)
    {
        # c(sum(degenesnum%in%u),length(u)) c(rownames(pwf)[u[-which(is.na(u))]])
        x <- u[which(u %in% DE_pwf)]
        x
    }, DE_pwf)
    
    xxx <- match2Genome(genome)
    
    pvals.6.gene.symbol <- sapply(pvals.6, function(u, xxx)
    {
        y <- xxx[match(u, as.character(xxx[, 2])), 1]
        y
    }, xxx)
    
    
    # Convert list to data frame
    pvals.6.df <- list_to_df(pvals.6)
    
    pvals.6.gene.symbol.df <- list_to_df(pvals.6.gene.symbol)
    
    dataset2 <- pvals.6.gene.symbol.df
    dataset2[sapply(dataset2, is.list)] <- sapply(dataset2[sapply(dataset2, 
        is.list)], function(x) sapply(x, function(y) paste(unlist(y), collapse = ", ")))
    
    temp.gene.name <- unique(apply(dataset2[, 2], 1, c))
    temp.gene.name.2 <- unique(gdata::trim(unlist(strsplit(temp.gene.name, split = ","))))
    
    DE_from_GO <- temp.gene.name.2
    
    colnames(pvals.6.df) <- c("category", "DEgene_ID")
    colnames(pvals.6.gene.symbol.df) <- c("category", "DEgene_symbol")
    
    # Finally, sort by p-value
    pvals <- pvals[order(pvals$over_represented_pvalue), ]
    
    # Supplement the table with the GO term name and ontology group but only if
    # the enrichment categories are actually GO terms
    if (any(grep("^GO:", pvals$category)))
    {
        GOnames <- select(GO.db, keys = pvals$category, columns = c("TERM", 
            "ONTOLOGY"))[, 2:3]
        colnames(GOnames) <- tolower(colnames(GOnames))
        pvals <- cbind(pvals, GOnames)
    }
    
    # And return
    pvals.2 <- merge(pvals, pvals.6.df, by = "category", sort = FALSE)
    
    pvals.3 <- merge(pvals.2, pvals.6.gene.symbol.df, by = "category", sort = FALSE)
    
    pvals.4 <- list(GO = pvals.3, DE_GO = DE_from_GO, cat2gene = cat2gene)
    
    return(pvals.4)
    
}

getGeneSet <- function(genes, genome, id, fetch.cats = c("GO:CC", "GO:BP", "GO:MF"))
{
    # Check for valid input
    if (any(!fetch.cats %in% c("GO:CC", "GO:BP", "GO:MF", "KEGG")))
    {
        stop("Invaled category specified.  Categories can only be GO:CC, GO:BP, GO:MF or KEGG")
    }
    # Convert from genome ID to org.__.__.db format
    orgstring <- as.character(.ORG_PACKAGES[match(gsub("[0-9]+", "", genome), 
        names(.ORG_PACKAGES))])
    # Multimatch or no match
    if (length(orgstring) != 1)
    {
        stop("Couldn't grab GO categories automatically.  Please manually specify.")
    }
    # Load the library
    library(paste(orgstring, "db", sep = "."), character.only = TRUE)
    # What is the default ID that the organism package uses?
    coreid <- strsplit(orgstring, "\\.")[[1]][3]
    
    # Now we need to convert it into the naming convention used by the organism
    # packages
    userid <- as.character(.ID_MAP[match(id, names(.ID_MAP))])
    # Multimatch or no match
    if (is.na(userid) | (length(userid) != 1))
    {
        stop("Couldn't grab GO categories automatically.  Please manually specify.")
    }
    # The (now loaded) organism package contains a mapping between the internal
    # ID and whatever the default is (usually eg), the rest of this function is
    # about changing that mapping to point from categories to the ID specified
    # Fetch the mapping in its current format Because GO is a directed graph, we
    # need to get not just the genes associated with each ID, but also those
    # associated with its children.  GO2ALLEGS does this.
    core2cat <- NULL
    if (length(grep("^GO", fetch.cats)) != 0)
    {
        # Get the name of the function which maps gene ids to go terms usually this
        # will be 'GO2ALLEG'
        gomapFunction <- .ORG_GOMAP_FUNCTION[orgstring]
        if (is.na(gomapFunction)) 
            gomapFunction <- .ORG_GOMAP_FUNCTION["default"]
        x <- toTable(get(paste(orgstring, gomapFunction, sep = "")))
        # Keep only those ones that we specified and keep only the names
        # core2cat=x[x$Ontology%in%gsub('^GO:','',fetch.cats),1:2]
        x[!x$Ontology %in% gsub("^GO:", "", fetch.cats), 2] <- "Other"
        core2cat <- x[, 1:2]
        colnames(core2cat) <- c("gene_id", "category")
    }
    if (length(grep("^KEGG", fetch.cats)) != 0)
    {
        x <- toTable(get(paste(orgstring, "PATH", sep = "")))
        # Either add it to existing table or create a new one
        colnames(x) <- c("gene_id", "category")
        if (!is.null(core2cat))
        {
            core2cat <- rbind(core2cat, x)
        } else
        {
            core2cat <- x
        }
    }
    
    # Now we MAY have to convert the 'gene_id' column to the format we are using
    if (coreid != userid)
    {
        # The mapping between user id and core id, don't use the <USER_ID>2<CORE_ID>
        # object as the naming is not always consistent
        user2core <- toTable(get(paste(orgstring, userid, sep = "")))
        # Throw away any user ID that doesn't appear in core2cat
        user2core <- user2core[user2core[, 1] %in% core2cat[, 1], ]
        # Make a list version of core2cat, we'll need it
        list_core2cat <- split(core2cat[, 2], core2cat[, 1])
        # Now we need to replicate the core IDs that need replicating
        list_core2cat <- list_core2cat[match(user2core[, 1], names(list_core2cat))]
        # Now we can replace the labels on this list with the user ones from
        # user2core, but there will be duplicates, so we have to unlist, label, then
        # relist
        user2cat <- split(unlist(list_core2cat, FALSE, FALSE), rep(user2core[, 
            2], sapply(list_core2cat, length)))
        # Now we only want each category listed once for each entry...
        user2cat <- sapply(user2cat, unique)
        ### In case you don't believe that this works as it should, here is the slow
        ### as all hell way for comparison... Make first list
        ### list_user2core=split(user2core[,1],user2core[,2]) Make the second
        ### list_core2cat=split(core2cat[,2],core2cat[,1]) Go through each entry in
        ### first list and expand using second...
        ### user2cat=sapply(list_user2core,function(u){unique(unlist(list_core2cat[u],FALSE,FALSE))})
        
    } else
    {
        # We don't need to convert anything (WOO!), so just make it into a list
        user2cat <- split(core2cat[, 2], core2cat[, 1])
        user2cat <- sapply(user2cat, unique)
    }
    # remove any empty strings
    user2cat <- lapply(user2cat, function(x)
    {
        if (length(x) > 1) 
            x = x[x != "Other"]
        x
    })
    
    ## we don't like case sensitivity
    names(user2cat) <- toupper(names(user2cat))
    gene2go <- user2cat[toupper(genes)]
    
    return(gene2go)
}

# Description: Prints progress through a loop copy from Matthew Young's
# goseq
pp <- function(total, count, i = i)
{
    if (missing(count))
    {
        count <- evalq(i, envir = parent.frame())
    }
    if (missing(total))
    {
        total <- evalq(stop, envir = parent.frame())
    }
    cat(round(100 * (count/total)), "%   \r")
}

outputGoBasedSelection <- function(Re.Go.adjusted.by.exon.SJ)
{
    
    # select GO term(10<=numInCat<=300 and BP only)
    
    index.select <- which(Re.Go.adjusted.by.exon.SJ[[1]]$numInCat >= 10 & Re.Go.adjusted.by.exon.SJ[[1]]$numInCat <= 
        300 & Re.Go.adjusted.by.exon.SJ[[1]]$ontology == "BP")
    
    Re.Go.adjusted.by.exon.SJ.select <- Re.Go.adjusted.by.exon.SJ[[1]][index.select, 
        ]
    Re.Go.adjusted.by.exon.SJ.select <- Re.Go.adjusted.by.exon.SJ.select[, -3]
    temp <- format(Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue, 
        scientific = TRUE, digits = 2)
    Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue <- temp
    
    rank.value.by.over_represented_pvalue <- rank(as.numeric(Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue), 
        ties.method = "min")
    
    Re.Go.adjusted.by.exon.SJ.select <- cbind(Re.Go.adjusted.by.exon.SJ.select, 
        rank.value.by.over_represented_pvalue)
    
    # Re.Go.adjusted.by.exon.SJ.select<-format(Re.Go.adjusted.by.exon.SJ.select,scientific
    # = TRUE,digits=2)
    
    return(Re.Go.adjusted.by.exon.SJ.select)
    
}

outputCatBasedSelection <- function(Re.Go.adjusted.by.exon.SJ)
{
    
    index.select <- which(Re.Go.adjusted.by.exon.SJ[[1]]$numInCat >= 10 & Re.Go.adjusted.by.exon.SJ[[1]]$numInCat <= 
        300)
    
    Re.Go.adjusted.by.exon.SJ.select <- Re.Go.adjusted.by.exon.SJ[[1]][index.select, 
        ]
    Re.Go.adjusted.by.exon.SJ.select <- Re.Go.adjusted.by.exon.SJ.select[, -3]
    temp <- format(Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue, 
        scientific = TRUE, digits = 2)
    Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue <- temp
    
    rank.value.by.over_represented_pvalue <- rank(as.numeric(Re.Go.adjusted.by.exon.SJ.select$over_represented_pvalue), 
        ties.method = "min")
    
    Re.Go.adjusted.by.exon.SJ.select <- cbind(Re.Go.adjusted.by.exon.SJ.select, 
        rank.value.by.over_represented_pvalue)
    
    # Re.Go.adjusted.by.exon.SJ.select<-format(Re.Go.adjusted.by.exon.SJ.select,scientific
    # = TRUE,digits=2)
    
    return(Re.Go.adjusted.by.exon.SJ.select)
    
}

# res11 <- getStaisitcs4Go(res1,gene.based.table)

getStaisitcs4Go <- function(GO.wall.DE_interest, mds.11.sample)
{
    
    GO.data <- GO.wall.DE_interest[[1]]
  
    y <- as.list(GO.data$DEgene_ID)
    
    re <- lapply(1:length(y), function(u, y, mds.11.sample)
    {
        
        yy <- y[[u]]
        
        y.id <- trim(c(unlist(strsplit(y[[u]], split = ","))))
        
        if (length(y.id) != 0)
        {
            
            yyy <- mean(as.numeric(unlist(mds.11.sample[match(y.id, mds.11.sample$geneID), 
                ]$numFeature)))
            
        } else
        {
            yyy <- 0
        }
        
        yyy
        
    }, y, mds.11.sample)
    
    re2 <- list_to_df(re)
    
    GO.data.1 <- cbind(GO.data, re2)
    GO.data.2 <- GO.data.1[, -(dim(GO.data.1)[2] - 1)]
    colnames(GO.data.2)[dim(GO.data.2)[2]] <- "Ave_value_DE"
    
    cat2gene <- GO.wall.DE_interest[[3]]
    
    rre <- lapply(1:length(cat2gene), function(u, cat2gene, mds.11.sample)
    {
        
        yy <- cat2gene[[u]]
        
        y.id <- yy
        
        if (length(y.id) != 0)
        {
            
            yyy <- mean(as.numeric(unlist(mds.11.sample[match(y.id, mds.11.sample$geneID), 
                ]$numFeature)), na.rm = TRUE)
            
        } else
        {
            yyy <- 0
        }
        
        yyy
        
    }, cat2gene, mds.11.sample)
    names(rre) <- names(cat2gene)
    rre2 <- list_to_df(rre)
    
    colnames(rre2) <- c("category", "Ave_value_all_gene")
    
    GO.data.3 <- merge(GO.data.2, rre2, by = "category", sort = FALSE)
    
    GO.data.3$Ave_value_DE <-  unlist(GO.data.3$Ave_value_DE)
    GO.data.3$Ave_value_all_gene <-  unlist(GO.data.3$Ave_value_all_gene)
    
    re3 <- list(GO.wall.DE_interest = GO.data.3, pwf.DE_interest = GO.wall.DE_interest[[2]])
    
    return(re3)
    
}

overlap_ratio <- function(x, y)
{
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x, y)))
}

color_scale <- function(c1 = "grey", c2 = "red")
{
    pal <- colorRampPalette(c(c1, c2))
    colors <- pal(100)
    return(colors)
}

getIdx <- function(v, MIN, MAX)
{
    if (MIN == MAX)
    {
        return(100)
    }
    intervals <- seq(MIN, MAX, length.out = 100)
    max(which(intervals <= v))
}

match2Genome <- function(genome_id)
{
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

reformatpath <- function(dir.name)
{
    CheckOPS <- Sys.info()[["sysname"]]
    
    if (CheckOPS == "Darwin")
    {
        temp <- unlist(strsplit(dir.name, split = "\\/"))
        
        if (!is.na(temp[3] == "H_driver"))
        {
            if (temp[3] == "H_driver")
            {
                temp[2] <- "Volumes"
                temp[3] <- "Bioinformatics$"
                dir.name <- do.call("file.path", as.list(temp))
            }
        }
        
    }
    
    return(dir.name)
}

getGeneSetBySize <- function(user2cat, go.size.limit)
{
    
    gene2go <- user2cat
    gene2go.select <- lapply(gene2go, function(x)
    {
        x = x[x != "Other"]
        x
    })
    
    gene2go.select.1 <- gene2go.select[lapply(gene2go.select, length) > 0]
    
    lower.size <- go.size.limit[1]
    upper.size <- go.size.limit[2]
    
    if (is.finite(lower.size) & is.finite(upper.size))
    {
        gene2go.select.2 <- gene2go.select.1[lapply(gene2go.select.1, length) > 
            lower.size & lapply(gene2go.select.1, length) <= upper.size]
    } else
    {
        
        gene2go.select.2 <- gene2go.select.1
    }
    
    return(gene2go.select.2)
    
}

# res <- PathwaySplice:::makeFeatureTable(res)
# 
makeFeatureTable <- function(jscs,use.multigene.aggregates = FALSE)
{
    temp <- fData(jscs)
    
    temp2 <- temp[which(temp$testable == TRUE), ]
    
    index.1 <- which(colnames(temp2) %in% c("geneID"))
    index.2 <- which(colnames(temp2) %in% c("countbinID"))
    
    index.3 <- which(colnames(temp2) %in% c("pvalue"))
    
    temp3 <- temp2[, c(index.1, index.2, index.3)]
    
    temp3 <- rapply(temp3, as.character, classes = "factor", how = "replace")
    
    if(use.multigene.aggregates == FALSE){
    temp3 <- temp3[-grep("\\+",temp3$geneID),]
    }
    
    row.names(temp3) = seq(1, dim(temp3)[1], 1)
    
    return(temp3)
}

nullpSplice = function(DEgenes, genome, id, bias.data = NULL, plot.fit = TRUE, 
    binsize = "auto")
    {
    # Input Checking
    if (!is.null(bias.data) & length(bias.data) != length(DEgenes))
    {
        stop("bias.data vector must have the same length as DEgenes vector!")
    }
    # Factors cause strange things to happen, remove them if they exist
    bias.data = unfactor(bias.data)
    DEgenes = unfactor(DEgenes)
    
    # Fetch length data from geneLenDataBase
    if (is.null(bias.data))
    {
        bias.data = getlength(names(DEgenes), genome, id)
    }
    
    # Fit a spline to the binary vector of DE calls vs gene length May not have
    # bias data for some of the entries, return NA at those positions
    pwf = rep(NA, length(DEgenes))
    w = !is.na(bias.data)
    pwf[w] = makespline(bias.data[w], DEgenes[w])
    
    # Make a data frame which contains all the data used to make the fit and the
    # fit itself
    out = data.frame(DEgenes = DEgenes, bias.data = bias.data, pwf = pwf, stringsAsFactors = FALSE)
    rownames(out) = names(DEgenes)
    
    # Plot the PWF if the arument has been specified
    if (plot.fit)
    {
        plotPwfSplice(out, binsize)
    }
    out
}

plotPwfSplice = function(pwf, binsize, pwf_col = 3, pwf_lwd = 2, xlab = "Biased Data in <binsize> gene bins.", 
    ylab = "Proportion of significant genes", ...)
    {
    # We shouldn't try and plot NAs obviously...
    w = !is.na(pwf$bias.data)
    o = order(pwf$bias.data[w])
    # What is the total range in the fit?
    rang = max(pwf$pwf, na.rm = TRUE) - min(pwf$pwf, na.rm = TRUE)
    if (rang == 0 & binsize == "auto") 
        binsize = 1000
    if (binsize == "auto")
    {
        # A low number of starting genes to bin, usually 100
        binsize = max(1, min(100, floor(sum(w) * 0.08)))
        resid = rang
        # Turn off warnings till we've worked out what we're doing
        oldwarn = options()$warn
        options(warn = -1)
        # Keep increasing the number of genes in each bin until the scatter around
        # the lines reaches the cutoff. Stop if we reach only 10 bins for the entire
        # plot
        while (binsize <= floor(sum(w) * 0.1) & resid/rang > 0.001)
        {
            binsize = binsize + 100
            # Assign each gene a 'bin number'
            splitter = ceiling(1:length(pwf$DEgenes[w][o])/binsize)
            # Determine the percentage DE in each bin
            de = sapply(split(pwf$DEgenes[w][o], splitter), mean)
            # Determine the average length in each bin
            binlen = sapply(split(as.numeric(pwf$bias.data[w][o]), splitter), 
                mean)
            # Calculate the residuals, how much the binned data deviates from the PWF
            resid = sum((de - approx(pwf$bias.data[w][o], pwf$pwf[w][o], binlen)$y)^2)/length(binlen)
        }
        options(warn = oldwarn)
    } else
    {
        # Assign each gene a 'bin number'
        splitter = ceiling(1:length(pwf$DEgenes[w][o])/binsize)
        # Determine the percentage DE in each bin
        de = sapply(split(pwf$DEgenes[w][o], splitter), mean)
        # Determine the average length in each bin
        binlen = sapply(split(as.numeric(pwf$bias.data[w][o]), splitter), mean)
    }
    # Now we've settled on a binsize, plot it Did the user specify the labels?
    # If so we can't put in the defaults or they'll be used twice and errors
    # result.
    xlab = gsub("<binsize>", as.character(binsize), xlab)
    if ("xlab" %in% names(list(...)))
    {
        if ("ylab" %in% names(list(...)))
        {
            plot(binlen, de, ...)
        } else
        {
            plot(binlen, de, ylab = ylab, ...)
        }
    } else if ("ylab" %in% names(list(...)))
    {
        plot(binlen, de, xlab = xlab, ...)
    } else
    {
        plot(binlen, de, xlab = xlab, ylab = ylab, ...)
    }
    # Add the PWF
    lines(pwf$bias.data[w][o], pwf$pwf[w][o], col = pwf_col, lwd = pwf_lwd)
}

# getResultsFromJunctionSeq
# 
# This function is used to get analysis results from using JunctionSeq
# 
# @param dir.name Path name for sample information file
# @param sample.file Sample information file
# @param count.file Count file
# @param gff.file Annotation file
# @param method.dispFinal Determine the method
# used to get a 'final' dispersion estimate.
# 
# @return The analysis result from JunctionSeq R package
# 
# @export
# 
# @examples
# 
# dir.name <- system.file('extdata', package='PathwaySplice')
# sample.file <- 'Sample_info.txt'
# count.file <- 'Counts.10.genes.txt'
# gff.file <- 'flat.chr22.10.genes.gff'
# res <- PathwaySplice:::getResultsFromJunctionSeq(dir.name, sample.file,
# count.file,gff.file, method.dispFinal = 'shrink',analysis.type = "exonsOnly")

getResultsFromJunctionSeq <- function(dir.name, sample.file, count.file, 
                                      gff.file, method.dispFinal = c("shrink", "max", "fitted", "noShare"),analysis.type)
{
    
            # set up method for calculating dispFinal
            method.dispFinal <- match.arg(method.dispFinal)
            
            # Get sample file
            dir.name <- reformatpath(dir.name)
              
            path.sample.file <- file.path(dir.name, sample.file)
            decoder.bySample <- read.table(path.sample.file, header = TRUE, stringsAsFactors = FALSE)
            
            x <- colnames(decoder.bySample)
  
            sample.ID.index <- which(colnames(decoder.bySample)==x[1])
            group.ID.index <- which(colnames(decoder.bySample)==x[2])
            Gender.index <- which(colnames(decoder.bySample)==x[3])
            
            # Get count file
            path.count.file <- file.path(dir.name, decoder.bySample[,sample.ID.index],count.file)
            
            # Get annotation file
            path.gff.file <- file.path(dir.name, "GTF_Files", gff.file)
                    
            # Analysis using exonsOnly,and adjust Gender
            jscs <- runJunctionSeqAnalyses(sample.files = path.count.file, sample.names = decoder.bySample[,sample.ID.index],condition = decoder.bySample[,group.ID.index], flat.gff.file = path.gff.file, analysis.type = analysis.type, nCores = 1, use.covars = decoder.bySample[, x[3], drop = FALSE], test.formula0 =formula(paste("~ ",paste("sample","countbin",paste0(x[3],":countbin"),sep="+"))), test.formula1 = formula(paste("~ ",paste("sample","countbin",paste0(x[3],":countbin"),"condition:countbin",sep = "+"))), effect.formula =formula(paste("~ ",paste("condition",x[3],"countbin",paste0(x[3],":countbin"),"condition:countbin",sep = "+")), geneLevel.formula = formula(paste("~ ",paste(x[3],"condition",sep = "+")), verbose = TRUE, debug.mode = TRUE, use.multigene.aggregates = TRUE,method.dispFinal = method.dispFinal)))
           
            return(jscs)
}

makeExample <- function(feature.table, num.gene)
{
  gene.name<- unique(feature.table$geneID)
  
  if(num.gene <= length(gene.name)){
  x <- sample(gene.name,num.gene)
  temp3 <- feature.table[which(feature.table$geneID %in% x),]
  row.names(temp3) <- seq(1,dim(temp3)[1])
  return(temp3)
  }else
  {
    cat("Please choose a number that is less or equal to the total number of genes\n")
  }
}

# res.reforamt <- PathwaySplice:::reformatPathwayOut(res)
reformatPathwayOut <- function(pathway.in){

  res <- dplyr::as_data_frame(pathway.in$GO)
  
  res
  
}

# PathwaySplice:::writeTibble(res.reforamt,"~/OutputTestPathwaySplice")
writeTibble <- function(tibble.input,output.dir){ 

  if (!dir.exists(output.dir))
  {
    dir.create(output.dir,recursive = TRUE)
  }
  
flatten_list = function(x){
  if (typeof(x) != "list") {
    return(x)
  }
  sapply(x, function(y) paste(y, collapse = " | "))
}

tibble.input %>%
  mutate_each(funs(flatten_list)) %>%
  write.csv(file.path(output.dir,"pathway.csv"))

}