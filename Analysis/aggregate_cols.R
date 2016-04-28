# GENERIC
aggregate_cols <- function(mat, col_groups = colnames(mat)) {
  agg_mat <- do.call(cbind, 
    lapply(
    	unique(col_groups), 
    	function(x) rowSums(as.matrix(mat[ , which(col_groups == x)]))
    	)
  )
  colnames(agg_mat) <- unique(col_groups)
  return(agg_mat)
}
