#' @noRd
discBinAreas <- function(object, object2, dist, resol, stype) {
  ad <- sqrt(object) / 2
  ad[2] <- sqrt(object2) / 2
  dAreas <- list()
  for (i in 1:2) {
    pt1 <- c(0, ifelse(i == 1, 0, dist))
    x1 <- pt1[1] - ad[i]
    x2 <- pt1[1] + ad[i]
    y1 <- pt1[2] - ad[i]
    y2 <- pt1[2] + ad[i]
    boun <- cbind(x = c(x1, x2, x2, x1, x1), y = c(y1, y1, y2, y2, y1))
    poly <- sf::st_sfc(sf::st_polygon(list(boun)))
    dAreas[[i]] <- sf::st_sample(
      poly,
      size = resol,
      type = stype,
      offset = c(0.5, 0.5)
    )
  }
  dAreas
}

#' @noRd
disc_bin_areas <- discBinAreas
