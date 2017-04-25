#' tiny.data
#' 
#' A tiny gene-based data set by sampling 500 genes from the whole gene list.
#' 
#' We used Hartley's QoRTs program to get counts of exons and/or splicing 
#'
#' junctions for each gene. The following is the count from genes 
#' 
#' 'ENSG00000135226+ENSG00000215110' for the sample SRR1660311:
#' 
#' ENSG00000135226+ENSG00000215110:A000
#' 
#' ENSG00000135226+ENSG00000215110:E001
#'
#' ENSG00000135226+ENSG00000215110:E002
#'
#' ENSG00000135226+ENSG00000215110:E003
#' 
#' ...
#' 
#' ENSG00000135226+ENSG00000215110:E005
#'
#' ENSG00000135226+ENSG00000215110:E006
#'
#' ENSG00000135226+ENSG00000215110:E007
#' 
#' ...
#' 
#' ENSG00000135226+ENSG00000215110:J016
#' 
#' ...
#' 
#' ESG00000135226+ENSG00000215110:J028
#'
#' We applied Hartley's JunctionSeq on counts to calculate Differential
#' 
#' Usage(DU) of exons and/or splicing junction, then converted 
#' 
#' DU results to gene based table. 
#' 
#' This table includes the following information:
#' 
#' \itemize{
#' 
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
#' In this gene based table, for 'ENSG00000135226+ENSG00000215110',
#'
#' we give row names as 'ENSG00000135226.ENSG00000135226+ENSG00000215110'
#'
#' and 'ENSG00000215110.ENSG00000135226+ENSG00000215110' to demonstrate 
#'
#' 'ENSG00000135226' and 'ENSG00000215110' are two genes from
#'
#' 'ENSG00000135226+ENSG00000215110'.
#'
#' @docType data
#' @keywords datasets
#' @name tiny.data
#' @usage data(TinyData)
#' 
#' @return A data set with 252 rows and 17 variables
#' 
#' @format A data frame with 252 rows and 17 variables
#' 
NULL