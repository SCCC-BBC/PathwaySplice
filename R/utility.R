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

reversemapping=function(map){
  tmp=unlist(map,use.names=FALSE)
  names(tmp)=rep(names(map),times=as.numeric(summary(map)[,1]))
  return(split(names(tmp),as.vector(tmp)))
}

