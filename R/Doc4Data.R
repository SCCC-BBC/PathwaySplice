#' PathwaySplice 
#'
#'@name PathwaySplice 
#'@doctype package 
#'@import  
NULL

#' hg19 Gene annotation data set
#' 
#' a data set for hg19 gene annotation
#'   
#' \itemize{
#'    \item hgnc_symbol. hgnc gene symbol  
#'    \item entrezgene. Entrez gene ID
#'    \item ensembl_gene_id. Ensembl gene ID
#'    \item chromosome_name. Chromosome information
#'    \item start_position. Starting position
#'    \item end_position. End position
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name hg19
#' @usage data(hg19)
#' @return a data frame with 64437 rows and 6 variables
#' @format A data frame with 64437 rows and 6 variables
NULL

#' An example data set
#' 
#' a data set for example data set
#' 
#' \itemize{
#'   \item geneID. Gene ID
#'   \item chr. Chromosome information
#'   \item start. Starting position
#'   \item end. End position
#'   \item strand. Strand information
#'   \item baseMean. BaseMean based on all samples
#'   \item geneWisePadj. perGeneQvalue from DEXSeq package 
#'   \item mostSigID. The feature that have smallest p-value within gene  
#'   \item mostSigPadjust. p-value for the most significant feature 
#'   \item numExons. Number of exons within gene
#'   \item numKnown. Number of splicing junction within gene
#'   \item numNovel. Number of novel splicing junction within gene
#'   \item exonsSig. Number of significant exons
#'   \item knownSig. Number of significant splicing junction
#'   \item novelSig. Number of significant novel splicing junction
#'   \item numFeatures. A list that includes number of exons,number of 
#'   splicing junction and number of novel splicing junction within gene 
#'   \item numSig. A list that includes number of significant exons,number of 
#'   significant splicing junction and number of significant novel splicing 
#'   junction within gene 
#' }
#'   
#' @docType data
#' @keywords datasets
#' @name mds
#' @usage data(mds)
#' @return a data set 23520 rows and 17 variables
#' @format A data frame with 23520 rows and 17 variables
#' 
NULL

#' hg38 Gene annotation data set
#' 
#' a data set for hg38 gene annotation
#'
#' \itemize{
#'    \item gene_id. Gene symbol  
#'    \item entrezgene. Entrez gene ID
#'    \item ensembl_gene_id. Ensembl gene ID
#'    \item chromosome_name. Chromosome information
#'    \item start_position. Starting position
#'    \item end_position. End position
#' }
#'
#' @docType data
#' @keywords datasets
#' @name hg38
#' @usage data(hg38)
#' @return a data set 24323 rows and 6 variables
#' @format A data frame with 24323 rows and 6 variables
NULL

#' mm10 Gene annotation data set
#' 
#' a data set for mm10 gene annotation
#'  
#' \itemize{
#'    \item gene_id. Gene symbol  
#'    \item entrezgene. Entrez gene ID
#'    \item ensembl_gene_id. Ensembl gene ID
#'    \item chromosome_name. Chromosome information
#'    \item start_position. Starting position
#'    \item end_position. End position
#' }
#'     
#' @docType data
#' @keywords datasets
#' @name mm10
#' @usage data(mm10)
#' @return a data set 23440 rows and 6 variables
#' @format A data frame with 23440 rows and 6 variables
NULL
