#' @export
#' @rdname gDist
gDist.rtop <- function(object, params = list(), ...) {
  params <- getRtopParams(object$params, newPar = params, ...)
  if (params$debug.level > 1) {
    debug.level <- params$debug.level
  } else {
    debug.level <- 0
  }
  if (!"dObs" %in% names(object)) {
    object$dObs <- rtopDisc(object$observations, params = params, ...)
  }
  if (!"dPred" %in% names(object) && "predictionLocations" %in% names(object)) {
    object$dPred <- rtopDisc(object$predictionLocations, params = params, ...)
  }
  object$gDistObs <- gDist(
    object$dObs,
    object$dObs,
    debug.level = debug.level,
    params = params,
    ...
  )
  if ("dPred" %in% names(object)) {
    object$gDistPredObs <- gDist(
      object$dObs,
      object$dPred,
      debug.level = debug.level,
      params = params,
      ...
    )
  }
  if ("dPred" %in% names(object)) {
    object$gDistPred <- gDist(
      object$dPred,
      object$dPred,
      diag = TRUE,
      debug.level = debug.level,
      params = params,
      ...
    )
  }
  object
}


#' @export
#' @noRd
gDist.SpatialPolygonsDataFrame <- function(object, ...) {
  utop_stop_legacy_sp(object)
}


#' @export
#' @noRd
gDist.SpatialPolygons <- function(object, ...) {
  utop_stop_legacy_sp(object)
}

#' @export
#' @noRd
gDist.STSDF <- function(object, ...) {
  utop_stop_legacy_sp(object, target = "stars")
}


#' @export
#' @rdname gDist
gDist.sf <- function(object, object2 = NULL, ...) {
  dObs <- rtopDisc(object, ...)
  gDistObs <- gDist(dObs, ...)
  if (!is.null(object2)) {
    dPred <- rtopDisc(object2, ...)
    gDistPredObs <- gDist(dObs, dPred, ...)
    gDistPred <- gDist(dPred, dPred, diag = TRUE, ...)
    list(
      gDistObs = gDistObs,
      gDistPred = gDistPred,
      gDistPredObs = gDistPredObs
    )
  } else {
    list(gDistObs = gDistObs)
  }
}

#' @export
#' @rdname gDist
gDist.stars <- function(object, object2 = NULL, ...) {
  object <- utop_stars_support(object)
  if (!is.null(object2)) {
    object2 <- utop_as_sf(object2)
  }
  gDist(object, object2, ...)
}


#' @export
#' @rdname gDist
gDist.list <- function(
  object,
  object2 = NULL,
  diag = FALSE,
  debug.level = 0,
  ...
) {
  variogramModel <- list(model = "Gho", params = 0)
  if (debug.level == 1) {
    print("Creating Ghos distances. This can take some time")
  }
  if (inherits(object[[1]], "sf") || inherits(object[[1]], "sfc")) {
    gDist <- varMat(
      object,
      object2,
      diag = diag,
      variogramModel = variogramModel,
      debug.level = debug.level,
      ...
    )
  } else {
    # These are the discretized points for hypotetical areas from binned variograms
    gDist <- data.frame(c1 = c(rep(0, length(object))), c2 = 0, cb = 0)
    lAreas <- lapply(object, FUN = function(aa) aa[[1]])
    gDist[, 1] <- mapply(
      vred,
      lAreas,
      MoreArgs = list(vredTyp = "ind", variogramModel = variogramModel)
    )
    lAreas <- lapply(object, FUN = function(aa) aa[[2]])
    gDist[, 2] <- mapply(
      vred,
      lAreas,
      MoreArgs = list(vredTyp = "ind", variogramModel = variogramModel)
    )
    gDist[, 3] <- mapply(
      vred,
      object,
      MoreArgs = list(vredTyp = "ind", variogramModel = variogramModel)
    )
  }
  return(as.matrix(gDist))
}
