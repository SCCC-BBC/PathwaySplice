
#' Title
#'
#' @param sample.name
#' @param condition
#' @param use.covars
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/Volumes/Bioinformatics$/2015/Nimer_Cheng/"
#' file.sample="decoder.bySample.txt"
#' path.file.sample<-paste0(dir.name,file.sample)
#' decoder.bySample<-read.table(path.file.sample,header=T,stringsAsFactors = F)
#' print(decoder.bySample)
#' Re<-GetDesign(decoder.bySample$sample.ID,decoder.bySample$group.ID)
#'
GetDesign <- function(sample.names,condition,use.covars=NULL) {

  if(! is.factor(condition)){
    condition <- factor(condition, levels = sort(unique(condition)))
  }

  design <- data.frame(condition = condition)
  if(! is.null(use.covars)){
    message(paste0("> rJSA: using covars:"," ",date()))
    if(class(use.covars) != "data.frame"){
      stop(paste0("FATAL ERROR: use.covars must be a data.frame! Instead it appears to be: ",class(use.covars)))
    }
    for(i in 1:ncol(use.covars)){
      if(! is.factor(use.covars[[i]]) ){
        use.covars[[i]] <- factor(use.covars[[i]], levels = sort(unique(use.covars[[i]]))  )
      }
    }

    design <- data.frame(cbind(design,use.covars))
    for(i in 1:length(names(use.covars))){
      message(paste0("      covar: ",names(use.covars)[i]))
      message(paste0(c("      ",paste0(use.covars[,i],collapse=", "))))
    }
    names(design) <- c("condition",names(use.covars))
  }
  row.names(design) <- sample.names

  print(design)

}
