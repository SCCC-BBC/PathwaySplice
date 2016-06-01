#' Title
#'
#' @param countfiles
#' @param countdata
#' @param samplenames
#' @param design
#' @param flat.gff.file
#' @param test.formula1
#' @param analysis.type
#' @param nCores
#' @param use.exons
#' @param use.junctions
#' @param use.known.junctions
#' @param use.novel.junctions
#' @param use.multigene.aggregates
#' @param gene.names
#' @param verbose
#' @param method.countVectors
#' @param noDESeqMatrix
#'
#' @return
#' @export
#'
#' @examples
#' dir.name="/Volumes/Bioinformatics$/2015/Nimer_Cheng/"
#' file.sample="decoder.bySample.txt"
#' file.count="_junction_seq_new_gtf_7/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
#' file.gff="Mus_musculus.GRCm38.83.JunctionSeq.flat.gff"

#' path.file.sample<-paste0(dir.name,file.sample)
#' decoder.bySample<-read.table(path.file.sample,header=T,stringsAsFactors = F)
#' print(decoder.bySample)

#' #Get count file
#' path.file.count<-paste0(dir.name,decoder.bySample$sample.ID,file.count)
#' countFiles<-paste0(path.file.count)
#' print(countFiles)

#' #Get annotation file
#' flat.file.gff<-paste0(dir.name,file.gff)
#' print(flat.file.gff)
#'
#' method.countVectors = c("geneLevelCounts","sumOfAllBinsForGene","sumOfAllBinsOfSameTypeForGene")
#' method.countVectors <- match.arg(method.countVectors)
#'
#' jscs = readJunctionSeqCounts(countfiles = as.character(sample.files),
#' samplenames = sample.names,
#' design = Re,
#' flat.gff.file = flat.gff.file,
#'  verbose = verbose,
#'  use.junctions = use.junctions,
#'  use.novel.junctions = use.novel.junctions,
#'  use.known.junctions = use.known.junctions,
#'  use.exons = use.exons,
#'  use.multigene.aggregates = use.multigene.aggregates,
#'  nCores = nCores,
#'  method.countVectors = method.countVectors,
#'  test.formula1 = test.formula1,
#'  gene.names = gene.names
#'  )
#'
#'
#'
readJunctionSeqCounts <- function(countfiles = NULL, countdata = NULL,
                                  samplenames,  design,
                                  flat.gff.file=NULL,
                                  test.formula1 = formula(~ sample + countbin + condition : countbin),
                                  analysis.type = c("junctionsAndExons","junctionsOnly","exonsOnly"),
                                  nCores = 1,
                                  use.exons = NULL, use.junctions = NULL,
                                  use.known.junctions = TRUE,
                                  use.novel.junctions = TRUE,
                                  use.multigene.aggregates = FALSE,
                                  gene.names = NULL,
                                  verbose = TRUE,
                                  method.countVectors = c("geneLevelCounts","sumOfAllBinsForGene","sumOfAllBinsOfSameTypeForGene"),
                                  noDESeqMatrix = FALSE)
{
  method.countVectors <- match.arg(method.countVectors)

  if(isTRUE(verbose)) {
    message("-> STARTING readJunctionSeqCounts (",date(),")")
  }


  analysis.type <- match.arg(analysis.type)
  if(is.null(use.junctions) && is.null(use.exons)){
    if(analysis.type == "junctionsAndExons"){
      use.junctions <- TRUE
      use.exons <- TRUE
    } else if(analysis.type == "junctionsOnly"){
      use.junctions <- TRUE
      use.exons <- FALSE
    } else if(analysis.type == "exonsOnly"){
      use.junctions <- FALSE
      use.exons <- TRUE
    }
  } else {
    if(is.null(use.junctions) || is.null(use.exons)){
      stop(paste0("Illegal syntax! If parameter use.junctions or use.exons are used, then BOTH must be set!\n use.junctions = '",use.junctions,"', use.exons = '",use.exons,"'"))
    }

    if(use.junctions && use.exons){
      analysis.type <- "junctionsAndExons"
    } else if(use.junctions && (! use.exons)){
      analysis.type <- "junctionsOnly"
    } else  if((! use.junctions) && use.exons){
      analysis.type <- "exonsOnly"
    } else {
      stop("Illegal syntax! Parameters use.exons and use.junctions cannot both be false!")
    }
  }
  if(isTRUE(verbose)){
    message("---> RJSC; (v",packageVersion("JunctionSeq"),")")
    message("---> RJSC: samplenames: ",paste0(samplenames, collapse=","))
    message("---> RJSC: flat.gff.file: ",flat.gff.file)
    message("---> RJSC: use.exons:",use.exons)
    message("---> RJSC: use.junctions:",use.junctions)
    message("---> RJSC: use.novel.junctions:",use.novel.junctions)
  }


  if((is.null(countfiles) && is.null(countdata))){
    stop("Fatal error: Either countfiles OR countdata must be set! Both are null!")
  }
  if(  (!is.null(countfiles)) && (!is.null(countdata))   ){
    stop("Fatal error: Either countfiles OR countdata must be set! Both are non-null!")
  }

  stopifnot( class(design) == "data.frame" )

  for(i in 1:ncol(design)){
    if( ! is.factor(design[[i]])){
      stop("ERROR: design must be a data.frame composed entirely of factors!")
    }
  }

  if(! is.null(countfiles)){
    lf <- lapply( countfiles, function(x)
      read.table( x, header=FALSE,stringsAsFactors=FALSE ) )
  } else {
    lf <- countdata
  }

  if( !all( sapply( lf[-1], function(x) all( x$V1 == lf[1]$V1 ) ) ) )
    stop( "Count files have differing gene ID column." )
  if(isTRUE(verbose)) message("---> File read complete.");

  dcounts <- sapply( lf, `[[`, "V2" )
  rownames(dcounts) <- lf[[1]][,1]
  dcounts <- dcounts[ substr(rownames(dcounts),1,1)!="_", ]

  bin.type <- sapply( rownames(dcounts),
                      function(x){
                        substr(strsplit(x, ":",fixed=TRUE)[[1]][2],1,1)
                      })
  raw.geneID <- sapply( rownames(dcounts),
                        function(x){
                          strsplit(x, ":",fixed=TRUE)[[1]][1]
                        })

  if(isTRUE(verbose)) message(paste0("---> Extracted counts. Found ",dim(dcounts)[1]," features so far."));

  geneCountTable <- dcounts[bin.type == "A",, drop=FALSE]
  rownames(geneCountTable) <- sapply(strsplit(rownames(geneCountTable), ":"),"[[",1)
  colnames(geneCountTable) <- as.character(samplenames)
  use.bins <- bin.type != "A"

  if(isTRUE(verbose)) message(paste0("---> Extracted gene-level counts. Found: ",dim(geneCountTable)[1], " genes and aggregate-genes."))
  if(isTRUE(verbose)) message(paste0("---> Removed gene features. Found: ",sum(use.bins), " features to be included so far."))

  if(isTRUE(! use.junctions) && isTRUE(! use.exons)){
    stop("FATAL ERROR: At least one of: use.junctions or use.exons must be set to TRUE. Otherwise you've got no data to test!")
  }

  if(isTRUE(! use.junctions)){
    use.bins <- use.bins & bin.type != "J" & bin.type != "N"
    if(isTRUE(verbose)) message(paste0("---> Removed splice junction features. Found: ",sum(use.bins), " features to be included so far."))
  }
  if(isTRUE(! use.novel.junctions)){
    use.bins <- use.bins & bin.type != "N"
    if(isTRUE(verbose)) message(paste0("---> Removed novel splice junction features. Found: ",sum(use.bins), " features to be included so far."))
  }
  if(isTRUE(! use.known.junctions)){
    use.bins <- use.bins & bin.type != "J"
    if(isTRUE(verbose)) message(paste0("---> Removed known splice junction features. Found: ",sum(use.bins), " features to be included so far."))
  }
  if(isTRUE(! use.exons)){
    use.bins <- use.bins & bin.type != "E"
    if(isTRUE(verbose)) message(paste0("---> Removed exon features. Found: ",sum(use.bins), " features to be included so far."))
  }

  is.multiGene.aggregate <- grepl("+", raw.geneID, fixed=TRUE)
  multiGene.aggregate.IDs <- unique(raw.geneID[is.multiGene.aggregate & use.bins])
  multiGene.aggregate.ct <- length(multiGene.aggregate.IDs)
  multiGene.aggregate.geneCt <- length(unlist( strsplit(multiGene.aggregate.IDs, "+",fixed=TRUE)))

  if(isTRUE(verbose)) message("---> Note: ",sum(is.multiGene.aggregate[use.bins])," counting bins from overlapping genes")
  if(isTRUE(verbose)) message("--->          There are ",multiGene.aggregate.ct,    " multigene aggregates.")
  if(isTRUE(verbose)) message("--->          There are ",multiGene.aggregate.geneCt," genes that are part of an aggregate.")
  if(isTRUE(! use.multigene.aggregates)){
    use.bins <- use.bins & (! is.multiGene.aggregate)
    if(isTRUE(verbose)) message(paste0("---> Removed multigene-aggregate features. Found: ",sum(use.bins), " features to be included so far."))
  }


  if(isTRUE(verbose)) message(paste0("---> Final feature count: ",sum(use.bins), " features to be included in the analysis."))
  dcounts <- dcounts[use.bins,, drop=FALSE]
  bin.type <- bin.type[use.bins]

  if(isTRUE(verbose)) message("---> Extracted feature counts.");

  colnames(dcounts) <- as.character(samplenames)
  splitted <- strsplit(rownames(dcounts), ":")
  exons <- sapply(splitted, "[[", 2)
  genesrle <- sapply( splitted, "[[", 1)

  if(isTRUE(verbose)) message("---> counts complete.");



  if(! is.null(flat.gff.file)){
    if(isTRUE(verbose)) message("-----> reading annotation...");
    anno.data <- readAnnotationData(flat.gff.file)
    if(isTRUE(verbose)) message("-----> formatting annotation...");
    #featureName     featureType     chrom   start   end     strand  gene_id part_number     transcripts

    exoninfo <- data.frame(chr = anno.data$chrom, start = anno.data$start, end = anno.data$end, strand = anno.data$strand)
    if(isTRUE(verbose)) message("-----> initial generation...");
    rownames(exoninfo) <- anno.data$featureName

    transcripts <- anno.data$transcripts
    transcripts <- gsub("\\+", ";", transcripts)
    names(transcripts) <- rownames(exoninfo)

    matching <- match(rownames(dcounts), rownames(exoninfo))
    if(any(is.na(matching))){
      stop("FATAL ERROR! Annotation file appears to be missing information! Are you sure you are using the correct flattened annotation file, created by prepare_annotation_with_splices.py?")
    }
    if(isTRUE(verbose)) message("-----> creating jscs...");
    jscs <- newJunctionSeqCountSet(countData=dcounts, geneCountData = geneCountTable, design=design, geneIDs=genesrle, countbinIDs=exons, featureIntervals=exoninfo[matching,], transcripts=transcripts[matching])
    jscs@annotationFile <- flat.gff.file
    jscs@flatGffData <- anno.data

    jscs@flatGffGeneData <- readGeneInfo(flat.gff.file)
    jscs <- mapGeneNames(jscs, gene.names)
  } else {
    if(isTRUE(verbose)) message("-> FINISHED readJunctionSeqCounts (",date(),")");
    message("Warning: flat gff annotation not set (via parameter flat.gff.file)! While technically optional, running without the annotation data may make interpretation of the data difficult. Much of the plotting functionality will not work!")
    warning("Warning: flat gff annotation not set (via parameter flat.gff.file)! While technically optional, running without the annotation data may make interpretation of the data difficult. Much of the plotting functionality will not work!")
    jscs <- newJunctionSeqCountSet(countData=dcounts, geneCountData = geneCountTable, design=design, geneIDs=genesrle, countbinIDs=exons);

  }
  attr(jscs,"AltMethods") <- c(attr(jscs,"AltMethods"), method.countVectors = method.countVectors)
  attr(jscs,"CallStack") <- list(deparse(match.call()))

  jscs@analysisType <- analysis.type
  featureChar <- substr(fData(jscs)$countbinID,1,1)
  fData(jscs)[["featureType"]] <- ifelse(featureChar == "E","exonic_part",ifelse(featureChar == "J", "splice_site", "novel_splice_site"))
  varMetadata( featureData(jscs) )[ "featureType", "labelDescription" ] <- "The type of feature (exonic_part,, splice_site, or novel_splice_site)."
  pData(jscs)$countfiles <- countfiles
  if(isTRUE(verbose)) message("-----> generating count vectors... (",date(),")")
  jscs@countVectors <- getAllJunctionSeqCountVectors(jscs, nCores = nCores, method.countVectors); #use.alternate.method = use.alternate.method)
  if(isTRUE(verbose)) message("-----> count vectors generated (",date(),")");

  if(! noDESeqMatrix){
    if(isTRUE(verbose)) message("-----> generating DESeqDataSet... (",date(),")")
    jscs <- makeDESeqDataSetFromJSCS(jscs, test.formula1 = test.formula1)
    if(isTRUE(verbose)) message("-----> DESeqDataSet generated (",date(),")")
    fData(jscs)[["allZero"]] <- (rowSums(counts(jscs)) == 0) |
      (rowSums(counts(jscs@DESeqDataSet)[, colData(jscs@DESeqDataSet)$countbin == "others"]) ==0)
    mcols(jscs@DESeqDataSet)$allZero <- fData(jscs)[["allZero"]]
  }

  return(jscs)
  if(isTRUE(verbose)) message("-> FINISHED readJunctionSeqCounts (",date(),")");
}
