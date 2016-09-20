dist2df <- function(x){
# convert distance matrix to long format data frame
  m  <- as.matrix(x)
  xy <- t(combn(colnames(m), 2))
  df <- data.frame(xy, dist = m[xy])
  return(df)
}
