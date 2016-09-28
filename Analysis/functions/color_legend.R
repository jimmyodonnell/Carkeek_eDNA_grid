color_legend <- function(col, lev){
# originally named 'legend.col'
    opar <- par

    n <- length(col)

    bx <- par("usr")

    box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
    bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
    box.cy <- c(bx[3], bx[3])
    box.sy <- (bx[4] - bx[3]) / n

    xx <- rep(box.cx, each = 2)

    par(xpd = TRUE)
    for(i in 1:n){

    yy <- c(box.cy[1] + (box.sy * (i - 1)),
    box.cy[1] + (box.sy * (i)),
    box.cy[1] + (box.sy * (i)),
    box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])

    }
    par(new = TRUE)
    plot(0, 0, type = "n",
    ylim = c(min(lev), max(lev)),
    yaxt = "n", ylab = "",
    xaxt = "n", xlab = "",
    frame.plot = FALSE)
    axis(side = 4, las = 2, tick = FALSE, line = .25)
    par <- opar
}

# myDat <- data.frame(
          # x = runif(100),
          # y = runif(100),
          # z = rnorm(100))
 
# myDat$Z <- cut(myDat$z, 100, label = FALSE)
 
# colr <- rev(heat.colors(100))
 
# plot(y ~ x, data = myDat,
     # col = colr[Z], pch = 20)
 
# color_legend(col = colr, lev = myDat$z)
