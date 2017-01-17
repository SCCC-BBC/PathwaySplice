#' makeGeneWiseTable
#'
#' This function is used to generate the gene-based results
#'
#' @param jscs output from JunctionSeq
#' @param gene.list gene list
#' @param FDR.threshold threshold to define whether exon or splicing junction is differential used
#' @param verbose whether to print runing process
#'
#' @return A gene based table
#' @export
#'
#' @examples
#'
#' #all.gene.list<-unique(as.character(fData(Re.example)$geneID))
#' 
#' #make a tiny example data set
#' 
#' #choosed.gene.list<-sample(all.gene.list,500)
#' 
#' #re.example.gene.based<-makeGeneWiseTable(Re.example,
#' #gene.list=choosed.gene.list)
#' #tiny.data<-re.example.gene.based
#' #save(tiny.data,file="./data/TinyData.RData")
#' 
#' #make an example data using all genes
#' 
#' #all.gene.data<-makeGeneWiseTable(Re.example,
#' #gene.list=all.gene.list)
#' #mds2<-all.gene.data
#' 
#' #save(mds2,file="./data/mds2.RData")
#' 
#'
makeGeneWiseTable <- function(jscs, gene.list, FDR.threshold = 0.05, 
    verbose = TRUE) {
    if (verbose) 
        message("   Compiling data table. ", date())
    
    mainTable <- data.frame(geneID = as.character(gene.list), 
        stringsAsFactors = FALSE)
    row.names(mainTable) <- gene.list
    mainTable <- AnnotatedDataFrame(mainTable)
    varMetadata(mainTable)["geneID", "labelDescription"] <- "Gene Unique Identifier"
    
    noGenes <- (length(gene.list) == 0)
    
    if (!noGenes) {
        geneAnno <- as.data.frame(t(sapply(gene.list, function(g) {
            geneRows <- which(fData(jscs)$geneID == g)
            c(as.character(fData(jscs)$chr[geneRows[1]]), as.numeric(min(fData(jscs)$start[geneRows])), 
                as.numeric(max(fData(jscs)$end[geneRows])), as.character(fData(jscs)$strand[geneRows[1]]))
        })))
        colnames(geneAnno) <- c("chr", "start", "end", "strand")
    } else {
        geneAnno <- data.frame(chr = character(), start = numeric(), 
            end = numeric(), strand = character())
    }
    
    mainTable$chr <- as.character(geneAnno$chr)
    mainTable$start <- geneAnno$start
    mainTable$end <- geneAnno$end
    mainTable$strand <- geneAnno$strand
    varMetadata(mainTable)[c("chr", "start", "end", "strand"), 
        "labelDescription"] <- c("Gene chromosome", "Gene start", 
        "Gene end", "Gene strand")
    
    # message('2')
    geneBaseMeans <- if (noGenes) {
        numeric()
    } else {
        rowMeans(jscs@geneCountData[match(gene.list, rownames(jscs@geneCountData)), 
            , drop = FALSE]/sizeFactors(jscs))
    }
    mainTable$baseMean <- if (noGenes) {
        character()
    } else {
        sprintf("%.1f", geneBaseMeans)
    }
    varMetadata(mainTable)["baseMean", "labelDescription"] <- "Gene BaseMean (simple normalized mean read or read-pair count per sample)"
    
    if (!is.null(fData(jscs)$geneWisePadj)) {
        mainTable$geneWisePadj <- if (noGenes) {
            numeric()
        } else {
            sapply(gene.list, function(g) {
                geneRows <- which(fData(jscs)$geneID == g)
                min(fData(jscs)$geneWisePadj[geneRows], na.rm = TRUE)
            })
        }
        varMetadata(mainTable)["geneWisePadj", "labelDescription"] <- "Gene-level adjusted p-value. P-value for the hypothesis that one or more features are DU."
    }
    
    mainTable$mostSigID <- if (noGenes) {
        character()
    } else {
        sapply(gene.list, function(g) {
            geneRows <- which(fData(jscs)$geneID == g)
            fData(jscs)$countbinID[geneRows[which.min(fData(jscs)$padjust[geneRows])]]
        })
    }
    varMetadata(mainTable)["mostSigID", "labelDescription"] <- "Feature ID of the most singificant feature."
    
    mainTable$mostSigPadjust <- if (noGenes) {
        numeric()
    } else {
        sapply(gene.list, function(g) {
            geneRows <- which(fData(jscs)$geneID == g)
            fData(jscs)$padjust[geneRows[which.min(fData(jscs)$padjust[geneRows])]]
        })
    }
    mainTable$mostSigPadjust <- sprintf("%.3g", mainTable$mostSigPadjust)
    varMetadata(mainTable)["mostSigPadjust", "labelDescription"] <- "Adjusted p-value of the most singificant feature."
    
    gene.row.list <- if (noGenes) {
        integer()
    } else {
        lapply(gene.list, function(g) {
            which(fData(jscs)$geneID == g)
        })
    }
    
    mainTable$numExons <- if (noGenes) {
        integer()
    } else {
        sapply(gene.row.list, function(geneRows) {
            sum(fData(jscs)$featureType[geneRows] == "exonic_part", 
                na.rm = TRUE)
        })
    }
    varMetadata(mainTable)["numExons", "labelDescription"] <- "Number of distinct exonic regions belonging to the gene."
    
    mainTable$numKnown <- if (noGenes) {
        integer()
    } else {
        sapply(gene.row.list, function(geneRows) {
            sum(fData(jscs)$featureType[geneRows] == "splice_site", 
                na.rm = TRUE)
        })
    }
    varMetadata(mainTable)["numKnown", "labelDescription"] <- "Number of distinct known splice sites belonging to the gene."
    mainTable$numNovel <- if (noGenes) {
        integer()
    } else {
        sapply(gene.row.list, function(geneRows) {
            sum(fData(jscs)$featureType[geneRows] == "novel_splice_site", 
                na.rm = TRUE)
        })
    }
    varMetadata(mainTable)["numNovel", "labelDescription"] <- "Number of distinct novel splice sites belonging to the gene."
    
    mainTable$exonsSig <- if (noGenes) {
        integer()
    } else {
        sapply(gene.row.list, function(geneRows) {
            sum(fData(jscs)$padjust[geneRows] < FDR.threshold & 
                fData(jscs)$featureType[geneRows] == "exonic_part", 
                na.rm = TRUE)
        })
    }
    varMetadata(mainTable)["exonsSig", "labelDescription"] <- paste0("Number of signficant exonic regions at p-adjust < ", 
        FDR.threshold)
    mainTable$knownSig <- if (noGenes) {
        integer()
    } else {
        sapply(gene.row.list, function(geneRows) {
            sum(fData(jscs)$padjust[geneRows] < FDR.threshold & 
                fData(jscs)$featureType[geneRows] == "splice_site", 
                na.rm = TRUE)
        })
    }
    varMetadata(mainTable)["knownSig", "labelDescription"] <- paste0("Number of signficant known splice junctions at p-adjust < ", 
        FDR.threshold)
    mainTable$novelSig <- if (noGenes) {
        integer()
    } else {
        sapply(gene.row.list, function(geneRows) {
            sum(fData(jscs)$padjust[geneRows] < FDR.threshold & 
                fData(jscs)$featureType[geneRows] == "novel_splice_site", 
                na.rm = TRUE)
        })
    }
    varMetadata(mainTable)["novelSig", "labelDescription"] <- paste0("Number of signficant novel splice junctions at p-adjust < ", 
        FDR.threshold)
    
    mainTable$numFeatures = if (noGenes) {
        character()
    } else {
        paste0(mainTable$numExons, "/", mainTable$numKnown, "/", 
            mainTable$numNovel)
    }
    varMetadata(mainTable)["numFeatures", "labelDescription"] <- "Number exonic regions / num known SJ / num novel SJ"
    mainTable$numSig = if (noGenes) {
        character()
    } else {
        paste0(mainTable$exonsSig, "/", mainTable$knownSig, "/", 
            mainTable$novelSig)
    }
    varMetadata(mainTable)["numSig", "labelDescription"] <- "Number sig exonic regions / num sig known SJ / num sig novel SJ"
    
    re <- ReformatData(mainTable)
    
    return(re)
}
