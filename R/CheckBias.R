#' Cbs
#'
#' Cbs is to check possible bias factor using the method in goseq
#'
#' @param re.gene.based Gene based table converted from DU results
#' @param adjust The possible bias factor to be used(gene length or number of exons)
#' @param sub_feature use exon("E") or splicing junction("J") to get gene table 
#' @param threshold cutoff for perGeneQvalue calculated using the method in DEXSeq
#' @param genomeID which genome(mm10 or hg19)
#' @param geneID which type of gene ID("ensGene" or "geneSymbol")
#'
#' @return A data frame that includes gene ID, status of differential gene
#' and probability weight function
#'
#' @export
#' @examples
#' 
#' data(mds)
#' Cbs(mds,adjust='E',sub_feature='E',threshold=0.05)
#'
#'
Cbs <- function(re.gene.based, adjust = "GL", sub_feature = NULL, 
    threshold, genomeID, geneID) {
    
    Data4Goterm <- re.gene.based
    
    if (is.null(sub_feature)) {
        Data4Goterm.sub_feature <- Data4Goterm
    } else {
        Data4Goterm.sub_feature <- Data4Goterm[grep(sub_feature, 
            Data4Goterm[, 8]), ]
    }
    
    if (sub_feature == "J") {
        Data4Goterm.sub_feature.geneID.NumOfJunctions <- Data4Goterm.sub_feature[, 
            c(1, 11)]
    } else {
        Data4Goterm.sub_feature.geneID.NumOfJunctions <- Data4Goterm.sub_feature[, 
            c(1, 10)]
    }
    
    Data4Goterm.sub_feature.Sig <- Data4Goterm.sub_feature[which(Data4Goterm.sub_feature[, 
        7] < threshold), ]
    
    # GO term analysis using GOSeq
    All.gene.id.based.on.sub_feature <- unique(Data4Goterm.sub_feature[, 
        1])
    
    All.gene.id.index <- rep(0, length(All.gene.id.based.on.sub_feature))
    names(All.gene.id.index) = All.gene.id.based.on.sub_feature
    
    All.genes.based.on.Sig.sub_feature <- unique(Data4Goterm.sub_feature.Sig[, 
        1])
    gene.DE_interest <- as.integer(which(All.gene.id.based.on.sub_feature %in% 
        All.genes.based.on.Sig.sub_feature))
    
    All.gene.id.index[gene.DE_interest] <- 1
    
    gene.with.matched.junction <- which(Data4Goterm.sub_feature.geneID.NumOfJunctions[, 
        1] %in% c(names(All.gene.id.index)))
    num.junction.4.matched.gene <- as.numeric(Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction, 
        2])
    
    All.gene.id.index.2 <- All.gene.id.index
    
    if (adjust == "GL") {
        pwf.DE_interest = nullp(All.gene.id.index.2, genomeID, 
            geneID, plot.fit = TRUE)
    } else {
        pwf.DE_interest = nullp(All.gene.id.index.2, genomeID, 
            geneID, bias.data = num.junction.4.matched.gene, 
            plot.fit = TRUE)
    }
    
    return(pwf.DE_interest)
    
}