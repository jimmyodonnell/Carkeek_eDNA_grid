nsamp <- 9

mat <- matrix(
  data = rnorm(nsamp*10), 
  nrow = nsamp
)

coords <- as.matrix(expand.grid(c(-1,0,1), c(-1,0,1)))

distmat <- as.matrix(dist(mat))

expand.grid(1:nsamp, 1:nsamp)
graph_data <- with(expand.grid(1:nsamp, 1:nsamp), data.frame(v1 = Var2, v2 = Var1))
graph_data$magnitude <- as.vector(distmat)

combn(1:9, 2)
?graph
duplicates <- which(apply(graph_data, 1, function(x)
  any(duplicated(x))))
graph_data <- graph_data[-duplicates,]

g <- graph.data.frame(graph_data, directed = FALSE)

plot.igraph(
  x = g, 
  layout = coords, 
  edge.curved = seq(-0.2, 0.2, length = ecount(g)), 
  edge.color = hsv(0,0,0,alpha = scale01(graph_data[,3])), 
  vertex.color = hsv(0.6,1,1), 
  vertex.frame.color = NA, 
  add = TRUE, 
  vertex.label = NA
)
points(coords, pch = 24)
plot(coords, pch = 24)
?igraph.plotting
?plot.igraph