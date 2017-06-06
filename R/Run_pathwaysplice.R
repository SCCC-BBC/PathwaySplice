#' makeGeneTable 
#' 
#' This function convert feature table into gene based table. 
#'  
#' @param feature.table An object loaded using \code{data(featureBasedData)}.
#' @param sig.threshold Threshold used for setting whether gene is differential or not based on smallest p-value of all features within a gene  
#'
#' @return A dataframe with the following arguments as cloumns for each gene 
#'   \itemize{
#'   \item geneID: Gene ID
#'   \item countbinID: feature ID in the feature differential usage analysis
#'   \item pvalue: unadjusted p value of each feature in the feature differential usage analysis
#' }
#' 
#' @export
#'
#' @examples
#' data(featureBasedData)
#' res <- makeGeneTable(featureBasedData)
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
#' that compares distributions of gene features for significant genes and non-significant genes. 
#' 
#' @param jscs.genewise.object A dataframe with genewise p-value for each gene, returned from \code{makeGeneTable()}
#' @param loc.x x coordinate for position of logistic regression p-value in figure
#' @param loc.y y coordinate for position of logistic regression p-value in figure
#' @param y.lim The largest number of exons in y axis in boxplot
#' @param boxplot.width Parameter for boxplot width
#'   
#' @details The logistic regression model Pr(significant gene) ~ number of features within the gene 
#' is implemented. Here features refer to exon bins or splicing junction bins, depending on 
#' how \code{genewise.pvalue} was obtained
#'   
#'       
#' @return NULL
#' @export
#' 
#' @examples
#' gene.based.table <- makeGeneTable(featureBasedData)
#' res <- lrTestBias(gene.based.table,loc.x=2,loc.y=150,y.lim=200,boxplot.width=0.3)
#' 
lrTestBias <- function(jscs.genewise.object, loc.x = 2, loc.y = 70, y.lim = 80, 
    boxplot.width = 0.3)
    {
    
    mydata <- jscs.genewise.object
    
    n.gene <- dim(mydata)[1]
    
    DE.out <- ifelse(mydata$sig.gene == 1, "Significant genes", "Non-significant genes")
    
    mydata.2 <- cbind(mydata, DE.out)
    
    par(mfrow = c(1, 1))
    
    if (var(as.numeric(unlist(mydata.2$numFeature))) != 0)
    {
        
        mylogit.2 <- glm(DE.out ~ as.numeric(numFeature), data = mydata.2, family = "binomial")
        re <- summary(mylogit.2)
        pvalue <- re$coefficients[2, 4]
        pvalue <- format(pvalue, width = 8, digits = 4)
        
        index.1 <- which(colnames(mydata.2) %in% c("numFeature"))
        index.2 <- which(colnames(mydata.2) %in% c("DE.out"))
        
        temp <- data.frame(mydata.2[, c(index.1, index.2)])
        
        temp$DE.out <- factor(temp$DE.out)
        
        temp$DE.out <- factor(temp$DE.out, levels = levels(temp$DE.out)[c(2, 
            1)])
        
        boxplot(unlist(temp$numFeature) ~ unlist(temp$DE.out), boxwex = boxplot.width, 
            ylab = "Number of features", col = "lightgray", ylim = c(1, y.lim))
        
        text(x = loc.x, y = loc.y, labels = c("", paste0("P-value from logistic regression:\n\n", 
            pvalue)), col = c(NA, "black"))
    } else
    {
        cat("There are no variations on the number of features\n")
    }
    
}

#' runPathwaySplice
#'
#' This function uses gene-based
#' table converted from makeGeneTable as an input, and select gene sets of certain type and size  
#' to perform gene set enrichment analysis adjusted by number of features in gene based table
#'  
#' @param res Gene based table 
#' @param genome Genome to be used(hg19 or mm10) 
#' @param id GeneID to be used(entrezgene or ensembl_gene_id)
#' @param gene2cat Get sets defined by users   
#' @param test.cats Gene set if users does not define their gene set 
#' @param go.size.cut Size of gene sets to be defined
#' @param method Method to be used for calculating overrepresented p value
#'        of gene sets(Options include Wallenius,Sampling, and Hypergeometric) 
#' @param repcnt Number of sampling if user use sampling method to calculate gene set enrichment p value
#' @param use.genes.without.cat Whether to include gene without mapping to get set for calculate gene set enrichment p value. if set FALSE, use the genes that have mapped gene sets only as background; if set TRUE, use all genes as background   
#'    
#'
#' @return A list that has gene set enrichment analysis results
#' @export
#'
#' @examples
#' res <- makeGeneTable(featureBasedData)
#' res <- runPathwaySplice(res,genome='hg19',id='ensGene',
#'                          test.cats=c('GO:BP'),
#'                          go.size.cut=c(5,30),
#'                          method='Wallenius')
#' 
#'  
runPathwaySplice <- function(res, genome, id, gene2cat = NULL, 
                             test.cats = c("GO:CC", "GO:BP", "GO:MF"),
                             go.size.cut = c(lower.size = 0, upper.size = NULL),
                             method = "Wallenius",
                             repcnt = 2000, use.genes.without.cat = FALSE)
{
  x <- res[, 3]
  names(x) <- res[, 1]
  pwf <- nullp(x, genome, id, bias.data = res[, 5], plot.fit = TRUE)
  CatDE <- pathwaysplice(pwf, genome = genome, id = id, gene2cat = gene2cat, 
                         test.cats = test.cats, go.size.cut = go.size.cut, method = method, repcnt = repcnt, 
                         use.genes.without.cat = use.genes.without.cat)
  CatDE
  
}

#' enrichmentMap
#'
#' enrichmentMap is used to draw Enrichment Map based on similarities defined using Jaccard Coefficient between GOs or gene sets
#'                                  
#' @param goseqres Object returned from runPathwaySplice
#' @param n Maximum number of category to be shown
#' @param fixed If set to FALSE, will invoke tkplot
#' @param vertex.label.font Font size of vertex label
#' @param similarity.threshold Threshold for defining Jaccard Coefficient(JC)
#'        
#'        JC ranges from 0 to 1:
#'        
#'        JC=0, indicates there are no overlap on genes between
#'               two gene sets
#'        
#'        JC=1, indicates two gene sets are identical  
#'        
#'        similarity.threshold=0, indicates the enrichment map includes
#'        all gene sets with their mutual JC greater than 0
#'        
#' @param output.file.dir Output dir for the gene set information file on network
#' @param label.vertex.by.index Which way to be used for labeling vertex on network
#'        
#'        FALSE indicates to label vertex by the name of gene sets
#'        
#'        TRUE indicates to label vertex by the index of gene sets    
#'          
#' @param ... Additional parameter 
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
#'                          go.size.cut=c(5,30),
#'                          method='Wallenius')
#' 
#' dir.name <- tempdir()
#' output.file.dir <- file.path(dir.name,'OutputEnmapEx')
#' 
#' enmap <- enrichmentMap(res,n=3,similarity.threshold=0,
#'                        output.file.dir = output.file.dir,
#'                        label.vertex.by.index = TRUE)
#'                        
enrichmentMap <- function(goseqres, n = 50, fixed = TRUE, 
                          vertex.label.font = 1, 
                          similarity.threshold,
                          output.file.dir, 
                          label.vertex.by.index = FALSE, ...)
    {
    
    if (!dir.exists(output.file.dir))
    {
        dir.create(output.file.dir)
    }
    
    GO.name <- goseqres[[1]]$category
    temp <- goseqres[[1]]$DEgene_ID
    names(temp) <- GO.name
    x <- goseqres[[1]]
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
        
        V(g)$name <- vertexname
        
        V(g)$color <- "red"
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
        E(g)$width <- sqrt(wd[, 3] * 5)
        
        g <- delete.edges(g, E(g)[wd[, 3] < similarity.threshold])
        
        idx <- unlist(sapply(V(g)$name, function(x) which(x == vertexname[1:n])))
        
        cols <- color_scale("red", "#E5C494")
        
        V(g)$color <- cols[sapply(pvalue, getIdx, min(pvalue), max(pvalue))]
        
        Edata <- as.data.frame(get.edgelist(g))
        Edata$edgewidth <- E(g)$width
        Vdata <- data.frame(pathway = V(g)$name, color = V(g)$color)
        map_data <- list(edge_data = Edata, vertex_data = Vdata)
        
        cnt <- as.integer(y$numDEInCat)
        
        names(cnt) <- vertexname[1:n]
        
        cnt2 <- cnt[V(g)$name]
        
        V(g)$size <- cnt2/sum(cnt2) * 100
    }
    
    netplot(g, vertex.label.font = vertex.label.font, vertex.label.color = "black", 
        fixed = fixed, ...)
    
    invisible(g)
    
    re2 <- map_data
    return(re2)
}

#' gmtGene2Cat
#'
#' Read a gene set file in GMT format, and return a list with its name
#' being a gene id, and each element of this list
#' being the pathways that this gene corresponds to
#'
#' @param dir.name Directory for the gene sets in GMT format that is located in 
#' @param pathway.file Input file for the gene sets in GMT format
#' @param file.type Indicates the gene set file in GMT format is in local or url
#' @param gene.anno.file Gene annotation file supplied as a file 
#' @param genomeID Genome ('mm10','hg19' or 'hg38') to be used
#'
#' @return A list with its names being geneID, its elements being the pathways
#'
#' @export
#'
#' @examples
#' 
#' dir.name <- system.file('extdata', package='PathwaySplice')
#' canonical.pathway.file <- '10.cp.gmt.txt'
#' res <- gmtGene2Cat(dir.name,canonical.pathway.file,'local',genomeID='hg19')
#' 
gmtGene2Cat <- function(dir.name, pathway.file, file.type,
                        gene.anno.file = NULL, 
                        genomeID = c("mm10", "hg19", "hg38"))
    {
    
    gmt_input_file <- file.path(dir.name, pathway.file)
    
    gene.2.cat.gmt <- gene2cat2(gmt_input_file, file.type)
    
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
#' res <- gmtGene2Cat(dir.name,canonical.pathway.file,
#'                    'local',genomeID='hg19')
#' gene.based.table <- makeGeneTable(featureBasedData)
#' 
#' res1 <- runPathwaySplice(gene.based.table,genome='hg19',
#'                          id='ensGene',gene2cat=res,
#'                          method='Wallenius')
#' 
#' res2 <- runPathwaySplice(gene.based.table,genome='hg19',
#'                          id='ensGene',gene2cat=res,
#'                          method='Hypergeometric')
#' 
#' dir.name <- tempdir()
#' output.dir <- file.path(dir.name,'OutputPostAnalysis')
#' 
#' output.file.name.1 <- 'In_ad_not_un.xls'
#' output.file.name.2 <- 'In_un_not_ad.xls'
#' res3 <- compareResults(4,res1,res2,output.dir,output.dir,
#'                       type.boxplot='Only3',
#'                       output.file.name.1,output.file.name.2)
#' @export
compareResults <- function(n.go, adjusted, unadjuasted, venn.dir, 
                          boxplot.dir,type.boxplot = c("All", "Only3"),
                          In.ad.not.un.file, In.un.not.ad.file)
{
  
  if (!dir.exists(venn.dir))
  {
    dir.create(venn.dir)
  }
  
  if (!dir.exists(boxplot.dir))
  {
    dir.create(boxplot.dir)
  }
  
  n <- n.go
  
  example.go.adjusted.by.exon <- adjusted
  example.go.unadjusted <- unadjuasted
  
  if (is.na(example.go.adjusted.by.exon) || is.na(example.go.unadjusted))
  {
    cat("One of results is empty\n\n")
    return()
  } else
  {
    
    if (dim(example.go.adjusted.by.exon$GO)[1] >= n && dim(example.go.unadjusted$GO)[1] >= 
        n)
    {
      
      adjusted <- example.go.adjusted.by.exon$GO[1:n, 1]
      unadjusted <- example.go.unadjusted$GO[1:n, 1]
      
      re <- list(adjusted = adjusted, unadjusted = unadjusted)
      
      venn.plot <- venn.diagram(x = re[c(1, 2)], filename = file.path(venn.dir, 
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
        index1 <- match(In.adjusted.not.in.unadjusted, example.go.adjusted.by.exon$GO$category)
        In.ad.not.un <- example.go.adjusted.by.exon$GO[index1, 
                                                                ]$Ave_value_all_gene
        
        yy <- cbind(example.go.unadjusted$GO[index1, ]$rank.value.by.over_represented_pvalue, 
                    example.go.adjusted.by.exon$GO[index1, ]$rank.value.by.over_represented_pvalue)
        
        
        index2 <- match(In.unadjusted.not.in.adjusted, example.go.unadjusted$GO$category)
        In.un.not.ad <- example.go.unadjusted$GO[index2, ]$Ave_value_all_gene
        
        yyy <- cbind(example.go.unadjusted$GO[index2, ]$rank.value.by.over_represented_pvalue, 
                     example.go.adjusted.by.exon$GO[index2, ]$rank.value.by.over_represented_pvalue)
        
        rre <- list(yy = yy, yyy = yyy)
        
        xx <- cbind(unlist(In.ad.not.un), unlist(In.un.not.ad))
        
        colnames(xx) <- c("In.ad.not.un", "In.un.not.ad")
        
        cp.top.adjusted.25 <- unlist(example.go.adjusted.by.exon$GO[1:n, 
                                                                             ]$Ave_value_all_gene)
        cp.top.unadjusted.25 <- unlist(example.go.unadjusted$GO[1:n, 
                                                                         ]$Ave_value_all_gene)
        
        cp.all.adjusted <- unlist(example.go.adjusted.by.exon$GO$Ave_value_all_gene)
        cp.all.unadjusted <- unlist(example.go.unadjusted$GO$Ave_value_all_gene)
        
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
          yy <- rbind(cbind(xx[, 1], rep("In.ad.not.un", length(xx[, 
                                                                   1]))), cbind(xx[, 2], rep("In.un.not.ad", length(xx[, 2]))), 
                      cbind(cp.top.adjusted.25, rep("cp.top.adjusted.25", length(cp.top.adjusted.25))), 
                      cbind(cp.top.unadjusted.25, rep("cp.top.unadjusted.25", 
                                                      length(cp.top.unadjusted.25))), cbind(cp.all.adjusted, 
                                                                                            rep("cp.all.adjusted", length(cp.all.adjusted))), cbind(cp.all.unadjusted, 
                                                                                                                                                    rep("cp.all.unadjusted", length(cp.all.unadjusted))))
          colnames(yy) <- c("y", "grp")
          yy <- as.data.frame(yy)
          png(file.path(boxplot.dir, "boxplot.png"))
          boxplot(as.numeric(as.character(y)) ~ grp, data = yy)
          dev.off()
        })
        
        Output_file <- file.path(boxplot.dir, In.ad.not.un.file)
        writegototable(example.go.adjusted.by.exon$GO[index1, 
                                                               ], Output_file)
        
        Output_file <- file.path(boxplot.dir, In.un.not.ad.file)
        writegototable(example.go.unadjusted$GO[index2, ], 
                       Output_file)
        
        return(rre)
        
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