#' featureBasedData
#' 
#' This dataset includes analysis results of RNA-seq data in Dolatshad et al. (2015), 
#' which compared transcriptome of CD34+ cells from myelodysplastic syndrome (MDS) patients 
#' with SF3B1 mutations vs. healthy controls using RNA sequencing. 
#' The \code{JunctionSeq} package was used to assess differential usage of counting bins, which are
#' non-overlapping segments of the exons or splicing junctions (see Fig 1 in Anders et al. (2012)). 
#' Because of the size limit, only counting bins associated with 
#' a subset of genes were included here for demonstration. 
#' 
#' @docType data
#' @keywords datasets
#' @name featureBasedData
#' @usage data(featureBasedData)
#'  
#' @format A data frame with variables for gene identifier (\code{geneID}), gene feature identifier (\code{countbinID}), and 
#'  p-value for gene feature (\code{pvalue}). Here we used "gene feature" and "counting bin" interchangeably
#' 
#' @references H Dolatshad, A Pellagatti, M Fernandez-Mercado1, B H Yip, L Malcovati, M Attwood, B Przychodzen
#' N Sahgal, A A Kanapin, H Lockstone, L Scifo, P Vandenberghe, E Papaemmanuil, C W J Smith, P J Campbell, 
#' S Ogawa1, J P Maciejewski, M Cazzola, K I Savage1 and J Boultwood1 (2015) \emph{Disruption of SF3B1 results in deregulated 
#' expression and splicing of key genes and pathways in myelodysplastic syndrome hematopoietic stem and progenitor 
#' cells.} Leukemia (2015) 29, 1092-1103
#' 
#' Anders S, Reyes A, Huber W (2012) \emph{Dececting differential usage of exons from RNA-seq data.} 
#' Genome Research 22(10): 2008-2017
#' 
NULL