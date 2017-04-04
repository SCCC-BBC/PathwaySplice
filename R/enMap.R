#' enrichmentmap
#'
#' enrichmentmap is used to draw network based on similarities between GOs
#'
#' @param GoSeqRes Object returned from Run_pathwaysplice
#' @param n Maximum number of category to be shown
#' @param fixed If set to FALSE, will invoke tkplot
#' @param vertex.label.font Font size of vertex label
#' @param SimilarityThreshold Threshold for defining Jaccard Coefficient(JC)
#'        
#'        JC ranges from 0 to 1:
#'        
#'        JC=0, indicates there are no overlap on genes between
#'               two gene sets
#'        
#'        JC=1, indicates two gene sets are identical  
#'        
#'        SimilarityThreshold=0, indicates the enrichment map includes
#'        all gene sets with their mutual JC greater than 0

#' @param output.file.dir Output dir for the gene set information file on network
#' @param label_vertex_by_index Which way to be used for labeling vertex on network
#'        
#'        FALSE indicates to label vertex by the name of gene sets
#'        
#'        TRUE indicates to label vertex by the index of gene sets    
#'          
#' @param ... Additional parameter 
#' @export
#' @return A figure for visualizing enrichment network
#' 
#' @author Aimin created this funciton based on enrichMap function in G Yu's DOSE R package
#' 
#' @examples
#'
#' res <- runpathwaysplice(tiny.data,adjust='exon_SJ',sub_feature='E',
#' 0.05,genomeID='hg19',geneID='ensGene', method='Wallenius')
#' 
#' output.file.dir <- "~/OutputEnmap"
#' 
#' enmap <- enrichmentmap(res,n=10,SimilarityThreshold=0,
#' output.file.dir = output.file.dir,label_vertex_by_index = TRUE)

enrichmentmap <-
  function(GoSeqRes,
           n = 50,
           fixed = TRUE,
           vertex.label.font = 1,
           SimilarityThreshold,output.file.dir,label_vertex_by_index = FALSE,
           ...) 
    {

     if (!dir.exists(output.file.dir)) {
      dir.create(output.file.dir)
     }
    
      GO.name <- GoSeqRes[[1]]$category
      temp <- GoSeqRes[[1]]$DEgene_ID
      names(temp) <- GO.name
      x <- GoSeqRes[[1]]
      geneSets <- temp

    y <- as.data.frame(x)
    
    if (any(grep("^GO:", y$category))) {
      VertexName <- paste0(y$term, ":", y$numDEInCat)
      
      if(label_vertex_by_index == TRUE) {
      VertexName.index <- seq(1,length(VertexName))
      
      output.text <- as.data.frame(cbind(VertexName.index,VertexName))[1:n,]
      VertexName <- VertexName.index
      colnames(output.text) <- c("index","name")
      write.table(output.text,
                  file = file.path(output.file.dir,"enrichmap_GO.xls"),
                  quote = FALSE,
                  col.names = TRUE,
                  row.names = FALSE,sep = "\t")
      }
    } else
    {
      VertexName <- paste0(y$category, ":", y$numDEInCat)
      
      if(label_vertex_by_index == TRUE) {
      VertexName.index <- seq(1,length(VertexName))
      
      output.text <- as.data.frame(cbind(VertexName.index,VertexName))[1:n,]
      VertexName <- VertexName.index
      colnames(output.text) <- c("index","name")
      
      write.table(output.text,
                  file = file.path(output.file.dir,"enrichmap_pathway.xls"),
                  quote = FALSE,
                  col.names = TRUE,
                  row.names = FALSE,sep = "\t")
      }
    }
    
    if (nrow(y) < n) {
      n <- nrow(y)
    }
    y <- y[1:n, ]
    
    if (n == 0) {
      stop("no enriched term found...")
    } else if (n == 1) {
      g <- graph.empty(0, directed = FALSE)
      g <- add_vertices(g, nv = 1)
      
      V(g)$name <- VertexName
      
      V(g)$color <- "red"
    } else {
      pvalue <- as.numeric(y$over_represented_pvalue)
      
      id <- y[, 1]
      geneSets <- geneSets[id]
      
      n <- nrow(y) #
      w <- matrix(NA, nrow = n, ncol = n)
      colnames(w) <- rownames(w) <- VertexName[1:n]
      
      for (i in 1:n) {
        for (j in i:n) {
          w[i, j] <- overlap_ratio(geneSets[id[i]], geneSets[id[j]])
        }
      }
      
      wd <- melt(w)
      wd <- wd[wd[, 1] != wd[, 2], ]
      wd <- wd[!is.na(wd[, 3]), ]
      g <- graph.data.frame(wd[, -3], directed = FALSE)
      E(g)$width <- sqrt(wd[, 3] * 5)
      
      g <- delete.edges(g, E(g)[wd[, 3] < SimilarityThreshold])
      
      idx <-
        unlist(sapply(V(g)$name, function(x)
          which(x == VertexName[1:n])))
      
      cols <- color_scale("red", "#E5C494")
      
      V(g)$color <-
        cols[sapply(pvalue, getIdx, min(pvalue), max(pvalue))]

      Edata <- as.data.frame(get.edgelist(g))
      Edata$edgewidth <- E(g)$width
      Vdata <- data.frame(pathway = V(g)$name, color = V(g)$color)
      map_data <- list(edge_data = Edata, vertex_data = Vdata)
      
      cnt <- as.integer(y$numDEInCat)

      names(cnt) <- VertexName[1:n]
      
      cnt2 <- cnt[V(g)$name]

      V(g)$size <- cnt2 / sum(cnt2) * 100
    }
    
    netplot(
      g,
      vertex.label.font = vertex.label.font,
      vertex.label.color = "black",
      fixed = fixed,
      ...
    )
 
    invisible(g)
    
    re2 <- map_data
    return(re2)
  }

overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y)) / length(unique(c(x, y)))
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