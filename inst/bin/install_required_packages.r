.cran_packages <- c("BiasedUrn","roxygen2","reshape2","igraph","gdata","VennDiagram","RColorBrewer","statmod","httr","plotrix","yaml","DBI","XML","RCurl","RSQLite","matrixStats","Hmisc","locfit","fastmatch","shiny")

.bioc_packages <- c("GenomeInfoDbData","GenomeInfoDb","AnnotationHub","JunctionSeq","goseq",
                    "SummarizedExperiment","Biobase","DOSE",
                    "org.Hs.eg.db","BiocGenerics","AnnotationDbi",
                    "GO.db","geneLenDataBase","interactiveDisplayBase","AnnotationDbi","DO.db")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)