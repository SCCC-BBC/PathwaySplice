#' Title
#'
#' @param list_for_df
#'
#' @return
#' @export
#'
#' @examples
#' test.goseq3.df<-list_to_df(test.goseq3)
#'
list_to_df <- function(list_for_df) {
  list_for_df <- as.list(list_for_df)

  nm <- names(list_for_df)
  if (is.null(nm))
    nm <- seq_along(list_for_df)

  df <- data.frame(name = nm, stringsAsFactors = FALSE)
  df$value <- unname(list_for_df)
  df
}
