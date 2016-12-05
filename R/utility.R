# Utility functions for PathwaySplice

list_to_df <- function(list_for_df) {
    list_for_df <- as.list(list_for_df)
    
    nm <- names(list_for_df)
    if (is.null(nm)) 
        nm <- seq_along(list_for_df)
    
    df <- data.frame(name = nm, stringsAsFactors = FALSE)
    df$value <- unname(list_for_df)
    df
}

reversemapping = function(map) {
    tmp = unlist(map, use.names = FALSE)
    names(tmp) = rep(names(map), times = as.numeric(summary(map)[, 
        1]))
    return(split(names(tmp), as.vector(tmp)))
}

ReformatData <- function(re.PJ.gene.based) {
    re <- pData(re.PJ.gene.based)
    
    no.re.testable.index <- which(as.character(re$mostSigID) == 
        "character(0)")
    re2 <- re[-no.re.testable.index, ]
    
    All.gene.id.based.on.sub_feature <- unique(unlist(strsplit(re2[, 
        1], "\\+")))
    All.gene.id.index <- rep(0, length(All.gene.id.based.on.sub_feature))
    names(All.gene.id.index) = All.gene.id.based.on.sub_feature
    
    
    reformat.gene.p <- do.call(rbind, sapply(All.gene.id.based.on.sub_feature, 
        function(u, re2) {
            x <- re2[grep(u, re2[, 1]), -1]
            x <- as.data.frame(t(x))
            # colnames(x)<-colnames(Data4Goterm) x
        }, re2))
    
    re3 <- as.data.frame(reformat.gene.p)
    re3 <- cbind(All.gene.id.based.on.sub_feature, re3)
    colnames(re3)[1] = "geneID"
    
    return(re3)
    
}

# @examples dir.name='/media/H_driver/2016/Yang/MACS/MACS/'
# reformatPath(dir.name)

reformatPath <- function(dir.name) {
    CheckOPS = Sys.info()[["sysname"]]
    
    if (CheckOPS == "Darwin") {
        temp = unlist(strsplit(dir.name, split = "\\/"))
        temp[2] = "Volumes"
        temp[3] = "Bioinformatics$"
        dir.name = paste0(paste0(temp, collapse = "/"), "/")
    }
    
    return(dir.name)
}
