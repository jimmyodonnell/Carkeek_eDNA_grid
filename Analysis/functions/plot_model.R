plot_model <- function(model_list_item, pred_vec, 
  line_color = "indianred", 
  line_type  = 3, 
  band_color = hsv(h = 1, s = 1, v = 0.1, alpha = 0.2)
  ){
  if(length(model_list_item$conf) > 0) {
    x_bounds <- c(x_pred, rev(x_pred))
    conf     <- model_list_item$conf
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
    lines(x = x_pred, y = model_list_item$conf[,"fit"], 
        col = line_color, lwd = 3, lty = line_type)
}
