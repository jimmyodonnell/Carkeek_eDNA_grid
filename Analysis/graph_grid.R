
the_edges <- data.frame(t(combn(attributes(mydist)$Labels, 2)))


the_edges$dis <- apply(the_edges, 1, FUN = function(x) as.matrix(mydist)[x[1],x[2]])

coords <- cbind(my_metadata[,colname_xcoord], log(my_metadata[,colname_ycoord] + 100))
rownames(coords) <- my_metadata[, colname_env_sample]
grid_graph <- graph.data.frame(the_edges, directed = FALSE)
plot.igraph(
  x = grid_graph, 
  layout = coords, 
  edge.curved = seq(-0.3, 0.3, length = ecount(grid_graph)), 
  edge.color = hsv(0,0,0, alpha = 0.2*scale01(the_edges[,3])), 
  edge.width = 3,
  vertex.color = NA, 
  vertex.frame.color = NA, 
  # add = TRUE, 
  vertex.label = NA
)
newcoords <- apply(coords, 2, scale11)
#  points(x = 0, y = -1, col = rainbow(5), cex = 2, pch = 19)
points(
		x = newcoords,
		col = mycolors[mypam$pamobject$clustering[my_metadata[,colname_env_sample]]],
		pch = 19, #as.character(mypam$pamobject$clustering[my_metadata$env_sample_name]-1)
		cex = 4
	)

help.search("lower triangle")