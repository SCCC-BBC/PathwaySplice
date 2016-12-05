#' Cbs
#'
#' Cbs is to check possible bias factor using the method in goseq
#'
#' @param re.gene.based Gene based table
#' @param ad The possible bias factor
#' @param sub_feature The possible bias factor
#' @param threshold threshold used
#' @param genomeID which genome
#' @param geneID which type of gene ID
#' @param gene_model Gene model
#' @param method which method
#'
#' @return A data frame that includes gene ID, status of differential gene
#' and probability weight function
#'
#' @export
#'
#' @examples
#'
#' data(mds)
#' Cbs(mds,ad='E',sub_feature='E',threshold=0.05)
#'
#'
Cbs <- function(re.gene.based, ad = "GL", sub_feature = NULL, 
    threshold, genomeID, geneID, gene_model, method) {
    
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
    print(length(All.gene.id.index))
    
    gene.with.matched.junction <- which(Data4Goterm.sub_feature.geneID.NumOfJunctions[, 
        1] %in% c(names(All.gene.id.index)))
    num.junction.4.matched.gene <- as.numeric(Data4Goterm.sub_feature.geneID.NumOfJunctions[gene.with.matched.junction, 
        2])
    
    All.gene.id.index.2 <- All.gene.id.index
    
    print(All.gene.id.index.2)
    
    if (ad == "GL") {
        pwf.DE_interest = nullp(All.gene.id.index.2, genomeID, 
            geneID, plot.fit = TRUE)
    } else {
        pwf.DE_interest = nullp(All.gene.id.index.2, genomeID, 
            geneID, bias.data = num.junction.4.matched.gene, 
            plot.fit = TRUE)
    }
    
    return(pwf.DE_interest)
    
}