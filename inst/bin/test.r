library("JunctionSeq");

message("-----")
sessionInfo();
message("-----")

#dir.name <- "T:\\DCEG\\Projects\\CCSS\\steve\\scrap\\pathwaySplice\\PathwaySplice\\inst\\extdata\\"
#dir.name <- system.file("extdata",package = "PathwaySplice")
dir.name <- "~/PathwaySplice/inst/extdata/"
sample.file <- "Sample_info.txt"
count.file <- "QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
gff.file <- "flat.chr22.gff"

path.sample.file <- file.path(dir.name, sample.file)
decoder.bySample <- read.table(path.sample.file, header = TRUE, stringsAsFactors = FALSE)
path.count.file <- file.path(dir.name, decoder.bySample$sample.ID, count.file)
path.gff.file <- file.path(dir.name, "GTF_Files", gff.file)



jscs <- runJunctionSeqAnalyses(sample.files = path.count.file, 
                               sample.names = decoder.bySample$sample.ID, condition = decoder.bySample$group.ID, 
                               flat.gff.file = path.gff.file, analysis.type = "exonsOnly", 
                               nCores = 1, 
                               use.covars = decoder.bySample[, "Gender", drop = FALSE], 
                               test.formula0 = ~sample + countbin + Gender:countbin, 
                               test.formula1 = ~sample + countbin + Gender:countbin + condition:countbin, 
                               effect.formula = ~condition + Gender + countbin + Gender:countbin + condition:countbin, 
                               geneLevel.formula = ~Gender + condition, 
                               verbose = TRUE, 
                               debug.mode = TRUE, use.multigene.aggregates = TRUE, method.dispFinal = "max")

message("-----")
sessionInfo();
message("-----")
