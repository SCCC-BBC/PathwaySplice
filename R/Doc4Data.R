#' PathwaySplice 
#'
#'@name PathwaySplice 
#'@doctype package
#' 
NULL

#' tiny.data
#' 
#' A tiny gene-based data set by sampling 500 genes from the whole gene list
#' 
#' \itemize{
#'   \item geneID: Gene ID
#'   \item chr: Chromosome information
#'   \item start: Starting position
#'   \item end: End position
#'   \item strand: Strand information
#'   \item baseMean: BaseMean based on all samples
#'   \item geneWisePadj: perGeneQvalue from DEXSeq package 
#'   \item mostSigID: The feature that have smallest p-value within gene  
#'   \item mostSigPadjust: p-value for the most significant feature 
#'   \item numExons: Number of exons within gene
#'   \item numKnown: Number of splicing junction within gene
#'   \item numNovel: Number of novel splicing junction within gene
#'   \item exonsSig: Number of significant exons
#'   \item knownSig: Number of significant splicing junction
#'   \item novelSig: Number of significant novel splicing junction
#'   \item numFeatures: A list that includes number of exons, number of 
#'   splicing junction and number of novel splicing junction within gene 
#'   \item numSig: A list that includes number of significant exons, number of 
#'   significant splicing junction and number of significant novel splicing 
#'   junction within gene 
#' }
#'   
#' @docType data
#' @keywords datasets
#' @name tiny.data
#' @usage data(TinyData)
#' @return A data set 252 rows and 17 variables
#' @format A data frame with 252 rows and 17 variables
#' 
NULL

#' mds.11.sample
#' 
#' A gene-based data set converted from the differential usage 
#' analysis based on 11 subjuects(excluding 309)
#' 
#' 
#' \itemize{
#'   \item geneID: Gene ID
#'   \item chr: Chromosome information
#'   \item start: Starting position
#'   \item end: End position
#'   \item strand: Strand information
#'   \item baseMean: BaseMean based on all samples
#'   \item geneWisePadj: perGeneQvalue from DEXSeq package 
#'   \item mostSigID: The feature that have smallest p-value within gene  
#'   \item mostSigPadjust: p-value for the most significant feature 
#'   \item numExons: Number of exons within gene
#'   \item numKnown: Number of splicing junction within gene
#'   \item numNovel: Number of novel splicing junction within gene
#'   \item exonsSig: Number of significant exons
#'   \item knownSig: Number of significant splicing junction
#'   \item novelSig: Number of significant novel splicing junction
#'   \item numFeatures: A list that includes number of exons,number of 
#'   splicing junction and number of novel splicing junction within gene 
#'   \item numSig: A list that includes number of significant exons,number of 
#'   significant splicing junction and number of significant novel splicing 
#'   junction within gene 
#' }
#'   
#' @docType data
#' @keywords datasets
#' @name mds.11.sample
#' @usage data(mds11)
#' @return A data set with 23520 rows and 17 variables
#' @format A data frame with 23520 rows and 17 variables
#' 
NULL