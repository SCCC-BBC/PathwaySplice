

#' Subset a data frame of class "prepGOglm"
#' 
#' This is an S3 generic function for \code{[} with the class \code{prepGOglm}.
#' @param x a data frame with class \code{prepGOglm}
#' @param ... for future use
#' @param drop logical value
#' @export
#' @method [ prepGOglm
#' 
#' @return A subset of \code{x}.
#' @rdname subset_data_frame
##' @author Gu Mi \email{mig@@stat.oregonstate.edu}, Yanming Di
##' \email{diy@@stat.oregonstate.edu}
#' 
`[.prepGOglm` <- function(x, ..., drop=TRUE) {
   structure(NextMethod(), class="prepGOglm")
}
