#' @noRd
plot_utop_variogram_cloud <- function(x, ...) {
  data <- utop_variogram_data(x)
  if ("ord" %in% names(data)) {
    data$np <- data$ord
  }
  class(data) <- "variogramCloud"
  plot(data, ...)
}

S7::method(plot, UtopVariogramCloud) <- plot_utop_variogram_cloud
