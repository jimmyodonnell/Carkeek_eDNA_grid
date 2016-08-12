scale11 <- function(x) {
  zeroone <- (x-min(x)) / (max(x)-min(x))
  return(zeroone * 2 - 1)
}
