#' @noRd
adfunc <- function(sampleVariogram, observations, amul) {
  if (is.null(sampleVariogram)) {
    area <- utop_area(observations)
    # alternative is variogram
  } else {
    area <- c(sampleVariogram$a1, sampleVariogram$a2)
  }
  amax <- max(area)
  amin <- min(area)
  areas <- axTicks(
    1,
    c(amin / 5, amax * 10, amul),
    usr = c(log10(amin / 5) - 1, log10(amax) + 1),
    log = TRUE,
    nintLog = Inf
  )
  areas <- areas[(min(which(areas > amin)) - 1):(max(which(areas < amax)) + 1)]
  areas
}


#' @noRd
dfunc <- function(sampleVariogram, observations, dmul) {
  if (is.null(sampleVariogram)) {
    dmax <- sqrt(bbArea(utop_bbox(observations))) / 2
    dmin <- min(stats::dist(utop_centroid_coordinates(observations)))
  } else if (variogram_is_cloud(sampleVariogram)) {
    dmax <- max(sampleVariogram$dist)
    dmin <- min(sampleVariogram$dist)
  } else {
    dmax <- max(sampleVariogram$dist)
    dmin <- min(sampleVariogram$dist[sampleVariogram$np > 2])
  }
  if (dmin < dmax / 1e8) {
    dmin <- dmax / 1e8
  }
  dists <- axTicks(
    1,
    c(dmin / 5, dmax * 10, dmul),
    usr = c(log10(dmin / 5) - 1, log10(dmax) + 1),
    log = TRUE,
    nintLog = Inf
  )
  dists[(min(which(dists > dmin)) - 1):(max(which(dists < dmax)) + 1)]
}
