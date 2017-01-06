#' enrichmentMap
#'
#' enrichmentMap is used to draw network based on similarities between GOs
#'
#' @param GoSeqRes object returned from Run_pathwaysplice
#' @param gene.set.type whether you are interested in GO, KEGG, or other pathways
#' @param n maximum number of category to shown
#' @param fixed if set to FALSE, will invoke tkplot
#' @param vertex.label.font font size of vertex label
#' @param SimilarityThreshold threshold for defining similarity between GOs
#' @param ... additional parameter
#' @export
#' @return A figure for visualizing enrichment network
#' 
#' @author Aimin created this funciton based on enrichMap function in G Yu's DOSE R package
#' 
#' @examples
#'
#' data(mds)
#' Example.Go.adjusted.by.exon<-Run_pathwaysplice(mds,ad="exon_SJ",
#' sub_feature="E",0.05,genomeID="hg19",geneID="ensGene",
#' gene_model=hg19,method="Sampling")
#' re.w.adjusted<-enrichmentMap(Example.Go.adjusted.by.exon,n=5,SimilarityThreshold=0)

#' Example.Go.unadjusted<-Run_pathwaysplice(mds,ad="exon_SJ",
#' sub_feature="E",0.05,genomeID="hg19",geneID="ensGene",
#' gene_model=hg19,method="Hypergeometric")
#' re.w.unadjusted<-enrichmentMap(Example.Go.unadjusted,n=5,SimilarityThreshold=0)
#'

enrichmentMap <-
  function(GoSeqRes,
           gene.set.type = "GO",
           n = 50,
           fixed = TRUE,
           vertex.label.font = 1,
           SimilarityThreshold,
           ...) {

    if (gene.set.type == "GO") {
      GO.name <- GoSeqRes[[1]]$category
      temp <- GoSeqRes[[1]]$DEgene_ID
      names(temp) <- GO.name
      x = GoSeqRes[[1]]
      geneSets = temp
    } else{
      GO.name <- GoSeqRes$category
      temp <- GoSeqRes$DEgene_ID
      names(temp) <- GO.name
      x = GoSeqRes[[1]]
      geneSets = temp
    }
    
    y <- as.data.frame(x)
    
    if (any(grep("^GO:", y$category))) {
      VertexName <- paste0(y$term, ":", y$numDEInCat)
    } else
    {
      VertexName <- paste0(y$category, ":", y$numDEInCat)
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
          w[i, j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
        }
      }
      
      wd <- melt(w)
      wd <- wd[wd[, 1] != wd[, 2], ]
      wd <- wd[!is.na(wd[, 3]), ]
      g <- graph.data.frame(wd[, -3], directed = FALSE)
      E(g)$width = sqrt(wd[, 3] * 20)
      
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