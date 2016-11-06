plot_model <- function(confidence_df, pred_vec,
  line_color = "indianred",
  line_type  = 3,
  band_color = hsv(h = 1, s = 1, v = 0.1, alpha = 0.2)
  ){
  if(length(confidence_df) > 0) {
    x_bounds <- c(pred_vec, rev(pred_vec))
    conf     <- confidence_df
    y_bounds <- c(conf[,"lwr"], rev(conf[,"upr"]))
    polygon(
      x = x_bounds,
      y = y_bounds,
      col = band_color,
      border = NA
    )
  } else {
    warning("Object you pointed to does not have confidence bounds")
  }
    # add fit line
    lines(x = pred_vec, y = confidence_df[,"fit"],
        col = line_color, lwd = 3, lty = line_type)
}
