# getResultsFromJunctionSeq This function is used to get analysis
# results from using JunctionSeq @param dir.name Path name for sample
# information file @param sample.file Sample information file @param
# count.file Count file @param gff.file Annotation file @param ...
# Additional parameters(Define them based on runJunctionSeqAnalyses
# function in JuntionSeq) @return The analysis result from JunctionSeq
# R package @export @examples dir.name <- system.file('extdata',
# package='PathwaySplice') sample.file <- 'Sample_info.txt' count.file
# <- 'Counts.10.genes.txt' gff.file <- 'flat.chr22.10.genes.gff' res
# <-PathwaySplice:::getResultsFromJunctionSeq(dir.name, sample.file,
# count.file,gff.file, method.dispFinal = 'shrink',analysis.type =
# 'exonsOnly')

getResultsFromJunctionSeq <- function(dir.name, sample.file, count.file, 
                                      gff.file, ...)
{
  
  # Get sample file
  dir.name <- reformatpath(dir.name)
  
  path.sample.file <- file.path(dir.name, sample.file)
  decoder.bySample <- read.table(path.sample.file, header = TRUE, stringsAsFactors = FALSE)
  
  x <- colnames(decoder.bySample)
  
  sample.ID.index <- which(colnames(decoder.bySample) == x[1])
  group.ID.index <- which(colnames(decoder.bySample) == x[2])
  Gender.index <- which(colnames(decoder.bySample) == x[3])
  
  # Get count file
  path.count.file <- file.path(dir.name, decoder.bySample[, sample.ID.index], 
                               count.file)
  
  # Get annotation file
  path.gff.file <- file.path(dir.name, "GTF_Files", gff.file)
  
  # Analysis using runJunctionSeqAnalyse, and adjust Gender
  jscs <- runJunctionSeqAnalyses(sample.files = path.count.file, sample.names = decoder.bySample[, 
                                                                                                 sample.ID.index], condition = decoder.bySample[, group.ID.index], 
                                 flat.gff.file = path.gff.file, use.covars = decoder.bySample[, 
                                                                                              x[3], drop = FALSE], test.formula0 = formula(paste("~ ", paste("sample", 
                                                                                                                                                             "countbin", paste0(x[3], ":countbin"), sep = "+"))), test.formula1 = formula(paste("~ ", 
                                                                                                                                                                                                                                                paste("sample", "countbin", paste0(x[3], ":countbin"), "condition:countbin", 
                                                                                                                                                                                                                                                      sep = "+"))), effect.formula = formula(paste("~ ", paste("condition", 
                                                                                                                                                                                                                                                                                                               x[3], "countbin", paste0(x[3], ":countbin"), "condition:countbin", 
                                                                                                                                                                                                                                                                                                               sep = "+"))), geneLevel.formula = formula(paste("~ ", paste(x[3], 
                                                                                                                                                                                                                                                                                                                                                                           "condition", sep = "+"))), ...)
  
  return(jscs)
}

makeExample <- function(feature.table, num.gene)
{
  gene.name <- unique(feature.table$geneID)
  
  if (num.gene <= length(gene.name))
  {
    x <- sample(gene.name, num.gene)
    temp3 <- feature.table[which(feature.table$geneID %in% x), ]
    row.names(temp3) <- seq(1, dim(temp3)[1])
    return(temp3)
  } else
  {
    cat("Please choose a number that is less or equal to the total number of genes\n")
  }
}

heatmap_wPCA = function(Data, g_level = NULL)
{
  
  Data.pca = prcomp(t(Data))
  hmcol <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
  
  if (is.null(g_level))
  {
    type_level = 1:ncol(Data)
    col_level = "black"
    
    with(data.frame(Data.pca$x), scatter3D(PC1, PC2, PC3, colvar = NULL, type = "h", 
                                           ticktype = "detailed", bty = "b2", cex = 1, xlab = "PC 1", ylab = "PC 2", 
                                           zlab = "PC 3", theta = 40, phi = 40, pch = type_level, col = col_level, 
                                           main = "Principal component analysis"))
    
    
    
    with(data.frame(Data.pca$x), text3D(x = PC1, y = PC2, z = PC3, colnames(Data), 
                                        col = "black", add = TRUE, colkey = FALSE, cex = 0.5))
  } else
  {
    type_level = 1:ncol(Data)
    TEMP = factor(g_level)
    uniq_label = levels(TEMP)
    levels(TEMP) = hmcol[ceiling(seq(length.out = length(levels(TEMP)), from = 1, 
                                     to = 256))]
    col_level = as.character(TEMP)
    uniq_col = levels(TEMP)
    
    Data.pca = prcomp(t(Data))
    with(data.frame(Data.pca$x), scatter3D(PC1, PC2, PC3, colvar = NULL, type = "h", 
                                           ticktype = "detailed", bty = "b2", cex = 1, xlab = "PC 1", ylab = "PC 2", 
                                           zlab = "PC 3", theta = 40, phi = 40, pch = type_level, col = col_level, 
                                           main = "Principal component analysis"))
    
    legend("topright", legend = uniq_label, pch = type_level, col = uniq_col, 
           cex = 1, inset = c(0.02))
    
    with(data.frame(Data.pca$x), text3D(x = PC1, y = PC2, z = PC3, colnames(Data), 
                                        col = "black", add = TRUE, colkey = FALSE, cex = 0.5))
  }
}

getCount4EachBamUsingJobArray <- function(input.bam.dir, input.bam.pattern, 
                                          gtffile.gtf, output.file.dir)
{
  
  index <- system("echo $LSB_JOBINDEX", intern = TRUE)
  
  u <- as.integer(index)
  
  if (!dir.exists(output.file.dir))
  {
    dir.create(output.file.dir, recursive = TRUE)
  }
  
  cmd.java.1 = "module load java/1.8.0_60"
  
  cmd.java.2 = "export _JAVA_OPTIONS=\"-Xmx5G\""
  
  cmd.java.3 = "java -jar $HOME/NGS_tools/QoRTs/QoRTs_1.1.8/QoRTs.jar QC --noGzipOutput --keepMultiMapped --stranded"
  
  cmd = paste(cmd.java.1, cmd.java.2, cmd.java.3, sep = ";")
  
  bam.list <- list.files(input.bam.dir, pattern = input.bam.pattern, 
                         all.files = TRUE, recursive = TRUE, full.names = TRUE)
  
  x <- lapply(bam.list, function(u)
  {
    
    sample.name <- basename(dirname(u))
    
    cmd1 = paste(cmd, u, gtffile.gtf, file.path(output.file.dir, sample.name), 
                 sep = " ")
    cmd1
  })
  
  cmd2 <- x[[u]]
  
  cat(cmd2, "\n\n")
  
  system(cmd2)
}

createBsubJobArrayRfun <- function(Rfun, job.name, wait.job.name)
{
  x <- useJobArrayOnPegasus("parallel", "72:00", 16, 25000, 8, job.name, 
                            wait.job.name)
  xx <- paste(x, paste0("\"R -e ", paste0("'", Rfun, "'"), "\""), sep = " ")
  xx
}

useJobArrayOnPegasus <- function(job.option = c("general", "parallel", 
                                                "bigmem"), Wall.time, cores, Memory, span.ptile, job.name, wait.job.name = NULL)
{
  
  job.option <- match.arg(job.option)
  
  job.name.array <- job.name
  
  switch(job.option, parallel = {
    cmd0 = paste(Wall.time, "-n", cores, "-q parallel -R 'rusage[mem=", 
                 Memory, "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu", 
                 sep = " ")
  }, bigmem = {
    cmd0 = paste(Wall.time, "-n", cores, "-q bigmem -R 'rusage[mem=", 
                 Memory, "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu", 
                 sep = " ")
  }, general = {
    cmd0 = paste(Wall.time, "-n", cores, "-q general -R 'rusage[mem=", 
                 Memory, "] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu", 
                 sep = " ")
  })
  
  if (!is.null(wait.job.name))
  {
    cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\")\"", " -P bbc -J \"", 
                  job.name, paste0("\" -o %J.", job.name.array, ".log "), paste0("-e %J.", 
                                                                                 job.name.array, ".err -W"))
  } else
  {
    cmd1 = paste0("bsub -P bbc -J \"", job.name, paste0("\" -o %J.", 
                                                        job.name.array, ".log "), paste0("-e %J.", job.name.array, 
                                                                                         ".err -W"))
  }
  
  cmd = paste(cmd1, cmd0, sep = " ")
  
  return(cmd)
}

# R -e
# 'library(PathwaySplice);PathwaySplice:::processBamFile('/projects/scratch/bbc/Project/Pengzhang_data2015/Alignment_len60','STAR_out.sorted.bam$','~/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf','/scratch/projects/bbc/aiminy_project/peng_junction')'

# R -e
# 'library(PathwaySplice);PathwaySplice:::processBamFile('/projects/scratch/bbc/Project/Pengzhang_data2015/Alignment_len60','STAR_out.sorted.bam$','~/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf','/scratch/projects/bbc/aiminy_project/peng_junction/count_strand_based')'

processBamFile <- function(input.bam.dir, input.bam.pattern, gtffile.gtf, 
                           output.file.dir)
{
  
  bam.list <- list.files(input.bam.dir, pattern = input.bam.pattern, 
                         all.files = TRUE, recursive = TRUE, full.names = TRUE)
  
  n <- length(bam.list)
  
  job.name <- paste0("Count[1-", n, "]")
  
  Rfun1 <- "library(PathwaySplice);re <- PathwaySplice:::getCount4EachBamUsingJobArray("
  input <- input.bam.dir
  input.bam.pattern <- input.bam.pattern
  processed.gene.gtf <- gtffile.gtf
  output <- output.file.dir
  Rfun2 <- ")"
  
  Rinput <- paste0("\\\"", input, "\\\",", "\\\"", input.bam.pattern, 
                   "\\\",", "\\\"", processed.gene.gtf, "\\\",", "\\\"", output, "\\\"")
  Rfun <- paste0(Rfun1, Rinput, Rfun2)
  
  counting <- createBsubJobArrayRfun(Rfun, job.name, wait.job.name = NULL)
  
  system(counting)
}

# PathwaySplice:::getResultsFromJunctionSeq2('~/Dropbox (BBSR)/BBSR
# Team
# Folder/Aimin_Yan/peng/count_strand_based','Sample_info.txt','QC.spliceJunctionAndExonCounts.forJunctionSeq.txt','mouse_str.gff','shrink','junctionsAndExons','~/Dropbox
# (BBSR)/BBSR Team
# Folder/Aimin_Yan/peng/count_strand_based/Output_jscs')

getResultsFromJunctionSeq2 <- function(dir.name, sample.file, count.file, 
                                       gff.file, method.dispFinal = c("shrink", "max", "fitted", "noShare"), 
                                       analysis.type, output.file.dir)
{
  
  if (!dir.exists(output.file.dir))
  {
    dir.create(output.file.dir, recursive = TRUE)
  }
  
  # set up method for calculating dispFinal
  method.dispFinal <- match.arg(method.dispFinal)
  
  # Get sample file
  dir.name <- reformatpath(dir.name)
  
  path.sample.file <- file.path(dir.name, sample.file)
  decoder.bySample <- read.table(path.sample.file, header = TRUE, stringsAsFactors = FALSE)
  
  x <- colnames(decoder.bySample)
  
  sample.ID.index <- which(colnames(decoder.bySample) == x[1])
  group.ID.index <- which(colnames(decoder.bySample) == x[2])
  
  # Get count file
  path.count.file <- file.path(dir.name, decoder.bySample[, sample.ID.index], 
                               count.file)
  
  # Get annotation file
  path.gff.file <- file.path(dir.name, "GTF_Files", gff.file)
  
  jscs <- runJunctionSeqAnalyses(sample.files = path.count.file, sample.names = decoder.bySample[, 
                                                                                                 sample.ID.index], condition = decoder.bySample[, group.ID.index], 
                                 flat.gff.file = path.gff.file, analysis.type = analysis.type, nCores = 1, 
                                 verbose = TRUE, debug.mode = TRUE, use.multigene.aggregates = TRUE, 
                                 method.dispFinal = method.dispFinal)
  
  save(jscs, file = file.path(output.file.dir, "jscs.RData"))
  
  # return(jscs)
}

# R -e
# 'library(PathwaySplice);PathwaySplice:::submitJob4Jscs('/scratch/projects/bbc/aiminy_project/peng_junction/count_strand_based','Sample_info.txt','QC.spliceJunctionAndExonCounts.forJunctionSeq.txt','Mus_musculus.GRCm38.83.processed.sorted_stranded.gff','shrink','junctionsAndExons','/scratch/projects/bbc/aiminy_project/peng_junction/count_strand_based/Output_jscs')'

submitJob4Jscs <- function(dir.name, sample.file, count.file, gff.file, 
                           method.dispFinal, analysis.type, output.file.dir)
{
  
  job.name <- "RunJscs"
  
  Rfun1 <- "library(PathwaySplice);re <- PathwaySplice:::getResultsFromJunctionSeq2("
  
  input <- dir.name
  sample.file <- sample.file
  count.file <- count.file
  gff.file <- gff.file
  method.dispFinal <- method.dispFinal
  analysis.type <- analysis.type
  output <- output.file.dir
  
  Rfun2 <- ")"
  
  Rinput <- paste0("\\\"", input, "\\\",", "\\\"", sample.file, "\\\",", 
                   "\\\"", count.file, "\\\",", "\\\"", gff.file, "\\\",", "\\\"", 
                   method.dispFinal, "\\\",", "\\\"", analysis.type, "\\\",", "\\\"", 
                   output, "\\\"")
  
  Rfun <- paste0(Rfun1, Rinput, Rfun2)
  
  jscs <- createBsubJobArrayRfun(Rfun, job.name, wait.job.name = NULL)
  
  system(jscs)
}

makeGffFile <- function(input.gtf.file, stranded = c("yes", "no"), out.gff.dir)
{
  
  cmd.java.1 = "module load java/1.8.0_60"
  cmd.java.2 = "export _JAVA_OPTIONS=\"-Xmx5G\""
  
  # input.gtf.file <-
  # '~/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf'
  
  input.gtf.file.name <- tools::file_path_sans_ext(basename(input.gtf.file))
  stranded <- match.arg(stranded)
  
  switch(stranded, yes = {
    cmd.java.3 = "java -jar $HOME/NGS_tools/QoRTs/QoRTs_1.1.8/QoRTs.jar makeFlatGff --stranded"
    cmd = paste(cmd.java.1, cmd.java.2, cmd.java.3, sep = ";")
    cmd1 <- paste(cmd, input.gtf.file, file.path(out.gff.dir, paste0(input.gtf.file.name, 
                                                                     "_stranded.gff"), sep = " "))
  }, no = {
    cmd.java.3 = "java -jar $HOME/NGS_tools/QoRTs/QoRTs_1.1.8/QoRTs.jar makeFlatGff"
    cmd = paste(cmd.java.1, cmd.java.2, cmd.java.3, sep = ";")
    cmd1 <- paste(cmd, input.gtf.file, file.path(out.gff.dir, paste0(input.gtf.file.name, 
                                                                     ".gff"), sep = " "))
  })
  
  system(cmd1)
  
}

# R -e
# 'library(PathwaySplice);PathwaySplice:::submitJob4makeGffFile('~/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf','yes','/scratch/projects/bbc/aiminy_project/peng_junction/count_strand_based/GTF_Files')'

submitJob4makeGffFile <- function(input.gtf.file, stranded, out.gff.dir)
{
  
  if (!dir.exists(out.gff.dir))
  {
    dir.create(out.gff.dir, recursive = TRUE)
  }
  
  job.name <- "makeGff"
  
  Rfun1 <- "library(PathwaySplice);re <- PathwaySplice:::makeGffFile("
  
  Rinput <- paste0("\\\"", input.gtf.file, "\\\",", "\\\"", stranded, 
                   "\\\",", "\\\"", out.gff.dir, "\\\"")
  Rfun2 <- ")"
  
  Rfun <- paste0(Rfun1, Rinput, Rfun2)
  
  cmd.gff <- createBsubJobArrayRfun(Rfun, job.name, wait.job.name = NULL)
  
  system(cmd.gff)
}

adjustBystatistics1 <- function(gene.based.table)
{
  
  z.value <- qchisq(gene.based.table$geneWisePvalue, 1, lower.tail = FALSE)
  
  m <- lm(z.value ~ gene.based.table$numFeature)
  
  z.value.adjusted = mean(z.value) + residuals(m)
  
  p.new <- pchisq(z.value.adjusted, df = 1, lower.tail = FALSE)
  
  gene.based.table$geneWisePvalue <- p.new
  gene.based.table
}

adjustBystatistics2 <- function(gene.based.table)
{
  
  z.value <- qchisq(gene.based.table$geneWisePvalue, 1, lower.tail = FALSE)
  
  m <- loess(z.value ~ gene.based.table$numFeature, span = 0.3)
  
  z.value.adjusted = mean(z.value) + residuals(m)
  
  p.new <- pchisq(z.value.adjusted, df = 1, lower.tail = FALSE)
  
  gene.based.table$geneWisePvalue <- p.new
  gene.based.table
}

adjustBystatistics3 <- function(gene.based.table, degree.poly)
{
  
  z.value <- qchisq(gene.based.table$geneWisePvalue, 1, lower.tail = FALSE)
  
  m <- lm(z.value ~ poly(gene.based.table$numFeature, degree.poly))
  
  z.value.adjusted = mean(z.value) + residuals(m)
  
  p.new <- pchisq(z.value.adjusted, df = 1, lower.tail = FALSE)
  
  gene.based.table$geneWisePvalue <- p.new
  gene.based.table
}

adjustBystatistics4 <- function(gene.based.table, nKnots = 6)
{
  
  nKnots <- round(nKnots)
  
  z.value <- qchisq(gene.based.table$geneWisePvalue, 1, lower.tail = FALSE)
  
  m <- gam(z.value ~ s(gene.based.table$numFeature, k = nKnots, bs = "cr"))
  
  z.value.adjusted = mean(z.value) + residuals(m)
  
  p.new <- pchisq(z.value.adjusted, df = 1, lower.tail = FALSE)
  
  gene.based.table$geneWisePvalue <- p.new
  gene.based.table
}

# input.gtf.file <- /media/H_driver/Aimin_project/GTF_Files/Homo_sapiens.GRCh38.84.gtf
# out.gff.dir <- 
# 
# proceessGtf4makeGffFile(input.gtf.file,out.gff.dir) 
#
proceessGtf4makeGffFile <- function(input.gtf.file,out.gff.dir,use.cluster=NULL)
{
  #/media/H_driver/Aimin_project/GTF_Files/Homo_sapiens.GRCh38.84.gtf 
  
  cmd.sh.1 = "tail -n +6" 
  cmd.sh.2 = '| awk -F "\t" "{OFS="\t"; $1 = "chr"$1; print}" | awk -F"\t" "{OFS="\t"; if($1=="chrMT") $1="chrM"; print}" | sort -k1,1 -k4,4n' 
  cmd.sh.3 = ">" 
    
  #Homo_sapiens.GRCh38.84.processed.sorted.22.gtf
  
  input.gtf.file.name <- tools::file_path_sans_ext(basename(input.gtf.file))
 
  cmd1 <- paste(cmd.sh.1, input.gtf.file,cmd.sh.2, cmd.sh.3,file.path(out.gff.dir, paste0(input.gtf.file.name, 
                                                                     "_processed.gtf"), sep = " "))
  print(cmd1)
  
  system(cmd1)
  
}