plot_points_nolab <- function(x,y){
  par(mar = c(2,3,1,1))
  plot(
    x = x,
    y = y,
    ylim = c(0,1),
    xlim = c(0, max(x)),
    xaxt = "n",
    pch = 21,
    cex = 1,
    col = hsv(h = 0, s = 1, v = 0, alpha = 0.5),
    bg  = rgb(0,0,0,alpha = 0.1 ), 
    xlab = "",
	ylab = "",
    axes = FALSE,
    las = 1
  )
  axis(side = 1, at = c(0, 1000, 2000, 3000, 4000), 
    labels = c(0, "", 2000, "", 4000), lwd = 0, lwd.ticks = 1)
  axis(side = 2, lwd = 0, lwd.ticks = 1, las = 1)
  box()
}
