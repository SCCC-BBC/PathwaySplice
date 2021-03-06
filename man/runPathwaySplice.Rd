% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Run_pathwaysplice.R
\name{runPathwaySplice}
\alias{runPathwaySplice}
\title{runPathwaySplice}
\usage{
runPathwaySplice(genewise.table, genome, id, gene2cat = NULL,
  test.cats = c("GO:CC", "GO:BP", "GO:MF"), go.size.limit = c(10, 200),
  method = "Wallenius", repcnt = 2000, use.genes.without.cat = FALSE,
  binsize = "auto", output.file = tempfile())
}
\arguments{
\item{genewise.table}{data frame returned from function \code{makeGeneTable}}

\item{genome}{Genome to be used, options are 'hg19' or 'mm10'}

\item{id}{GeneID, options are 'entrezgene' or 'ensembl_gene_id'}

\item{gene2cat}{Get sets to be tested, these are defined by users, can be obtained from \code{gmtGene2Cat} function}

\item{test.cats}{Default gene ontology gene sets to be tested if \code{gene2cat} is not defined}

\item{go.size.limit}{Size limit of the gene sets to be tested}

\item{method}{the method used to calculate pathway enrichment p value. 
Options are 'Wallenius', 'Sampling', and 'Hypergeometric'}

\item{repcnt}{Number of random samples to be calculated when 'Sampling' is used, this argument
ignored unless \code{method='Sampling'}}

\item{use.genes.without.cat}{Whether genes not mapped to any gene_set tested are included in the analysis.
Default is set to FALSE, where genes not mapped to any tested categories are ignored in analysis. 
Set this option to TRUE if it's desired that all genes in \code{genewise.table} to be counted towards  
the  total  number  of  genes outside  the  category.}

\item{binsize}{The number of genes in each gene bin in the bias plot}

\item{output.file}{File name for the analysis result in .csv format.}
}
\value{
runPathwaySplice returns a \href{https://cran.r-project.org/package=dplyr}{tibble} with the following information:
\item{gene_set}{Name of the gene set. Note in this document we used the terms gene_set, category, 
and pathway interchangeably} 
\item{over_represented_pvalue}{P-vaue for the associated gene_set being over-represented among significant genes} 
\item{under_represented_pvalue}{P-vaue for the associated gene_set being under-represented among significant genes} 
\item{numSIGInCat}{The number of significant genes in the gene_set} 
\item{numInCat}{The total number of genes in the gene_set}                                          
\item{description}{Description of the gene gene_set} 
\item{ontology}{The domain of the gene ontology terms if GO categories were tested. 
      Go categories can be classified into three domains: cellular component, biological process, molecular function.} 
\item{SIGgene_ensembl}{Ensembl gene ID of significant genes in the gene_set}
\item{SIGgene_symbol}{Gene symbols of signficant genes in the gene_set}
\item{Ave_value_all_gene}{The average value for \code{numFeature} for all the genes in the gene_set, 
     note that \code{numFeature} is the bias factor adjusted by PathwaySplice}

These information are also saved in the file \code{output.file}
}
\description{
This function identifies pathways that are enriched with signficant genes, while accounting for 
different number of gene features (e.g. exons) associated with each gene
}
\details{
This function implements the methodology described in Young et al. (2011) to adjust for 
         different number of gene features (column \code{numFeature} in \code{gene.based.table}). 
         For example, gene features can be non-overlapping exon counting bins associated with each gene (Fig 1 in Anders et al. 2012).
         In the bias plot, the genes are grouped by \code{numFeature} in \code{genewise.table} into gene bins, 
         the proportions of signficant genes are then plotted against the gene bins.
}
\examples{
gene.based.table <- makeGeneTable(featureBasedData)

res <- runPathwaySplice(gene.based.table,genome='hg19',id='ensGene',
                         test.cats=c('GO:BP'),
                         go.size.limit=c(5,30),
                         method='Wallenius',binsize=20)
\dontrun{

# demonstrate how output file can be specified                 
res <- runPathwaySplice(gene.based.table,genome='hg19',id='ensGene',
                       test.cats=c('GO:BP'),
                       go.size.limit=c(5,30),
                       method='Wallenius',binsize=800, 
                       output.file=tempfile())    

# demonstrate using customized gene sets
dir.name <- system.file('extdata', package='PathwaySplice')
hallmark.local.pathways <- file.path(dir.name,'h.all.v6.0.symbols.gmt.txt')
hlp <- gmtGene2Cat(hallmark.local.pathways, genomeID='hg19')

res <- runPathwaySplice(gene.based.table,genome='hg19',id='ensGene',
                       gene2cat=hlp,
                       go.size.limit=c(5,200),
                       method='Wallenius',binsize=20, 
                       output.file=tempfile())
                       
}

}
\references{
Young MD, Wakefield MJ, Smyth GK, Oshlack A (2011) \emph{Gene ontology analysis for RNA-seq: 
accounting for selection bias}. Genome Biology 11:R14

Anders S, Reyes A, Huber W (2012) \emph{Dececting differential usage of exons from RNA-seq data.} 
Genome Research 22(10): 2008-2017
}
