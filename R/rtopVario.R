# Replace the dependent variable by its OLS residuals against the trend
# basis functions of formulaString, so that the sample variogram describes
# the residual (detrended) field used by the universal kriging system.
#' @noRd
ukDetrendForVariogram <- function(object, formulaString, params, discPoints) {
  if (!inherits(params, "rtopParams")) {
    params <- getRtopParams(params)
  }
  object$ukResidual <- ukResiduals(formulaString, object, params, discPoints)
  list(object = object, formulaString = ukResidual ~ 1)
}


#' @export
#' @rdname rtopVariogram
rtopVariogram.sf <- function(
  object,
  formulaString,
  params = list(),
  cloud,
  discPoints = NULL,
  ...
) {
  if (missing(object)) {
    stop("rtopVariogram: Observations are missing")
  }
  if (!inherits(params, "rtopParams")) {
    params <- getRtopParams(params, ...)
  }
  object <- utop_add_area(object)
  if (missing(formulaString)) {
    formulaString <- utop_default_formula(object)
    warning(paste("formulaString missing, using", formulaString))
  }
  if (!inherits(formulaString, "formula")) {
    formulaString <- as.formula(formulaString)
  }
  if (hasUkTrend(formulaString)) {
    detr <- ukDetrendForVariogram(object, formulaString, params, discPoints)
    object <- detr$object
    formulaString <- detr$formulaString
  }

  observations <- suppressWarnings(sf::st_centroid(object))
  observations$area <- utop_area(object)
  clvar <- gstat::variogram(formulaString, observations, cloud = TRUE, ...)
  .BigInt <- attr(clvar, ".BigInt")
  clvar$ord <- clvar$np
  clvar <- as.data.frame(clvar)
  clvar$a1 <- observations$area[clvar$left]
  clvar$a2 <- observations$area[clvar$right]

  if (missing(cloud)) {
    cloud <- params$cloud
  }
  if (cloud) {
    clvar$acl1 <- clvar$left
    clvar$acl2 <- clvar$right
    clvar <- clvar[, -which(names(clvar) %in% c("left", "right"))]
    var3d <- clvar
    class(var3d) <- c("rtopVariogramCloud", "data.frame")
    attr(var3d, ".BigInt") <- .BigInt
    var3d$np <- 1
  } else {
    abins <- adfunc(NULL, observations, params$amul)
    dbins <- dfunc(NULL, observations, params$dmul)
    observations$acl <- findInterval(observations$area, abins)
    clvar$acl1 <- observations$acl[clvar$left]
    clvar$acl2 <- observations$acl[clvar$right]

    ich <- which(clvar$acl1 > clvar$acl2)
    acl1c <- clvar$acl1
    clvar$acl1[ich] <- clvar$acl2[ich]
    clvar$acl2[ich] <- acl1c[ich]

    clvar$dbin <- findInterval(clvar$dist, dbins)
    clvar$np <- 1
    varnp <- aggregate(
      list(np = clvar$np),
      list(acl1 = clvar$acl1, acl2 = clvar$acl2, dbin = clvar$dbin),
      sum
    )
    var3d <- data.frame(
      np = varnp$np,
      aggregate(
        list(
          dist = clvar$dist,
          gamma = clvar$gamma,
          a1 = clvar$a1,
          a2 = clvar$a2
        ),
        list(acl1 = clvar$acl1, acl2 = clvar$acl2, dbin = clvar$dbin),
        mean
      )
    )
    class(var3d) <- c("rtopVariogram", "data.frame")
  }
  var3d
}


#' @export
#' @noRd
rtopVariogram.SpatialPolygonsDataFrame <- function(object, ...) {
  utop_stop_legacy_sp(object)
}


#' @export
#' @noRd
rtopVariogram.SpatialPointsDataFrame <- function(object, ...) {
  utop_stop_legacy_sp(object)
}

#' @export
#' @noRd
rtopVariogram.STSDF <- function(object, ...) {
  utop_stop_legacy_sp(object, target = "stars")
}


# Alternative binning:
#  x <- matrix(rnorm(30000), ncol=3)
#  breaks <- seq(-1, 1, length=5)
#  xints <- data.frame(
#  x1=cut(x[, 1], breaks=breaks),
#  x2=cut(x[, 2], breaks=breaks),
#  x3=cut(x[, 3], breaks=breaks))
#  table(complete.cases(xints))
#  xtabs(~ ., xints)

###############################
#' @export
#' @rdname rtopVariogram
rtopVariogram.rtop <- function(object, params = list(), ...) {
  params <- getRtopParams(object$params, newPar = params, ...)
  observations <- object$observations
  formulaString <- object$formulaString

  # calling the class-specific rtopVariogram method
  var3d <- rtopVariogram(
    observations,
    formulaString,
    params,
    discPoints = object$dObs,
    ...
  )
  if (inherits(var3d, "rtopVariogramCloud")) {
    object$variogramCloud <- var3d
  } else {
    object$variogram <- var3d
  }
  object
}


################################

#' @export
#' @rdname rtopVariogram
rtopVariogram.stars <- function(
  object,
  formulaString,
  params = list(),
  cloud,
  abins,
  dbins,
  discPoints = NULL,
  ...
) {
  if (!inherits(params, "rtopParams")) {
    params <- getRtopParams(params, ...)
  }
  amul <- params$amul
  dmul <- params$dmul
  if (missing(cloud)) {
    cloud <- params$cloud
  }
  if (missing(object)) {
    stop("rtopVariogram: Observations are missing")
  }
  observations <- utop_stars_add_area(object)
  support <- utop_stars_support(observations)
  if (missing(formulaString)) {
    formulaString <- utop_default_formula(observations)
    warning(paste("formulaString missing, using", formulaString))
  }
  if (!inherits(formulaString, "formula")) {
    formulaString <- as.formula(formulaString)
  }
  depvar <- as.character(formulaString[[2]])
  obs_mat <- utop_stars_attr_matrix(observations, depvar)

  if (hasUkTrend(formulaString)) {
    nspace <- utop_stars_nspace(observations)
    ntime <- utop_stars_ntime(observations)
    Fsp <- ukTrendMatrix(formulaString, support, params, discPoints)
    obs_mat <- matrix(
      ukResiduals(
        formulaString,
        depValues = as.vector(obs_mat),
        trendMatrix = Fsp[rep(seq_len(nspace), times = ntime), , drop = FALSE]
      ),
      nrow = nspace
    )
  }

  nspace <- utop_stars_nspace(observations)
  ntime <- utop_stars_ntime(observations)
  vmat <- matrix(0, nrow = nspace, ncol = nspace)
  indmat <- vmat

  for (ind in seq_len(ntime)) {
    vals <- obs_mat[, ind]
    findx <- which(!is.na(vals))
    nspace1 <- length(findx)
    if (nspace1 < 2) {
      next
    }
    vals <- vals[findx]
    ff <- outer(vals, vals, FUN = function(x, y) (y - x)^2 / 2)
    vmat[findx, findx] <- vmat[findx, findx] + ff
    indmat[findx, findx] <- indmat[findx, findx] + 1
  }
  vmat <- vmat / indmat
  dmat <- utop_dist_matrix(utop_centroid_coordinates(support))

  vario <- matrix(NA, ncol = 7, nrow = nspace * (nspace - 1) / 2)
  icount <- 0
  for (istat in 1:(nspace - 1)) {
    njs <- nspace - istat
    vario[(icount + 1):(icount + njs), ] <- matrix(
      c(
        dmat[istat, (istat + 1):nspace],
        vmat[istat, (istat + 1):nspace],
        rep(support$area[istat], njs),
        support$area[(istat + 1):nspace],
        rep(istat, njs),
        (istat + 1):nspace,
        indmat[istat, (istat + 1):nspace]
      ),
      ncol = 7
    )
    icount <- icount + njs
  }
  vario <- data.frame(vario)
  names(vario) <- c("dist", "gamma", "a1", "a2", "acl1", "acl2", "np")
  vario <- vario[!is.na(vario$gamma), ]

  if (cloud) {
    var3d <- vario
    class(var3d) <- c("rtopVariogramCloud", "data.frame")
  } else {
    if (missing(abins)) {
      abins <- adfunc(NULL, support, amul)
    }
    if (missing(dbins)) {
      dbins <- dfunc(NULL, support, dmul)
    }
    support$acl <- findInterval(support$area, abins)
    vario$acl1 <- support$acl[vario$acl1]
    vario$acl2 <- support$acl[vario$acl2]

    ich <- which(vario$acl1 > vario$acl2)
    acl1c <- vario$acl1
    vario$acl1[ich] <- vario$acl2[ich]
    vario$acl2[ich] <- acl1c[ich]

    vario$dbin <- findInterval(vario$dist, dbins)
    varnp <- aggregate(
      list(np = vario$np),
      list(acl1 = vario$acl1, acl2 = vario$acl2, dbin = vario$dbin),
      sum
    )
    var3d <- data.frame(
      np = varnp$np,
      aggregate(
        list(
          dist = vario$dist,
          gamma = vario$gamma,
          a1 = vario$a1,
          a2 = vario$a2
        ),
        list(acl1 = vario$acl1, acl2 = vario$acl2, dbin = vario$dbin),
        mean
      )
    )
    class(var3d) <- c("rtopVariogram", "data.frame")
  }
  var3d
}
