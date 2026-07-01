#' @noRd
uk_detrend_for_variogram <- function(
  object,
  formula,
  params,
  disc_points = NULL
) {
  object$ukResidual <- ukResiduals(formula, object, params, disc_points)
  list(object = object, formula = ukResidual ~ 1)
}

#' @noRd
bin_variogram_cloud <- function(clvar, observations, params) {
  abins <- adfunc(NULL, observations, params@amul)
  dbins <- dfunc(NULL, observations, params@dmul)
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
  data.frame(
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
}

#' @noRd
compute_variogram_sf <- function(
  object,
  formula,
  params,
  cloud = params@cloud,
  disc_points = NULL,
  ...
) {
  object <- utop_add_area(object)
  if (is.null(formula)) {
    formula <- as.formula(utop_default_formula(object))
    warning("formula missing, using ", deparse(formula), call. = FALSE)
  }
  if (!inherits(formula, "formula")) {
    formula <- as.formula(formula)
  }
  if (hasUkTrend(formula)) {
    detr <- uk_detrend_for_variogram(object, formula, params, disc_points)
    object <- detr$object
    formula <- detr$formula
  }

  observations <- suppressWarnings(sf::st_centroid(object))
  observations$area <- utop_area(object)
  clvar <- gstat::variogram(formula, observations, cloud = TRUE, ...)
  .BigInt <- attr(clvar, ".BigInt")
  clvar$ord <- clvar$np
  clvar <- as.data.frame(clvar)
  clvar$a1 <- observations$area[clvar$left]
  clvar$a2 <- observations$area[clvar$right]

  if (cloud) {
    clvar$acl1 <- clvar$left
    clvar$acl2 <- clvar$right
    clvar <- clvar[, -which(names(clvar) %in% c("left", "right"))]
    var3d <- clvar
    attr(var3d, ".BigInt") <- .BigInt
    var3d$np <- 1
    return(utop_wrap_variogram_result(var3d, cloud = TRUE))
  }

  utop_wrap_variogram_result(
    bin_variogram_cloud(clvar, observations, params),
    cloud = FALSE
  )
}

#' @noRd
compute_variogram_stars <- function(
  object,
  formula,
  params,
  cloud = params@cloud,
  abins = NULL,
  dbins = NULL,
  disc_points = NULL,
  ...
) {
  amul <- params@amul
  dmul <- params@dmul
  observations <- utop_stars_add_area(object)
  support <- utop_stars_support(observations)
  if (missing(formula) || is.null(formula)) {
    formula <- as.formula(utop_default_formula(observations))
    warning("formula missing, using ", deparse(formula), call. = FALSE)
  }
  if (!inherits(formula, "formula")) {
    formula <- as.formula(formula)
  }
  depvar <- as.character(formula[[2]])
  obs_mat <- utop_stars_attr_matrix(observations, depvar)

  if (hasUkTrend(formula)) {
    nspace <- utop_stars_nspace(observations)
    ntime <- utop_stars_ntime(observations)
    fsp <- ukTrendMatrix(formula, support, params, disc_points)
    obs_mat <- matrix(
      ukResiduals(
        formula,
        depValues = as.vector(obs_mat),
        trendMatrix = fsp[rep(seq_len(nspace), times = ntime), , drop = FALSE]
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
    return(utop_wrap_variogram_result(vario, cloud = TRUE))
  }

  if (is.null(abins)) {
    abins <- adfunc(NULL, support, amul)
  }
  if (is.null(dbins)) {
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
  utop_wrap_variogram_result(var3d, cloud = FALSE)
}

#' Create a sample variogram
#'
#' @param object A [Utop] object, `sf`, or `stars` data.
#' @param ... Method-specific arguments. For [Utop] objects, `params` updates
#'   parameters. For `sf` and `stars`, supply `formula`, `params`, and
#'   optional `disc_points`.
#'
#' @return A [Utop] object with a variogram attached, or a standalone
#'   [UtopVariogram] / [UtopVariogramCloud] object.
#' @export
utop_variogram <- S7::new_generic(
  name = "utop_variogram",
  dispatch_args = "object",
  fun = function(object, ...) {
    S7::S7_dispatch()
  }
)

S7::method(utop_variogram, Utop) <- function(object, params = list(), ...) {
  object@params <- utop_params(
    params = object@params,
    new_params = params,
    observations = object@observations,
    formula = object@formula
  )

  vario <- utop_variogram(
    object@observations,
    formula = object@formula,
    params = object@params,
    disc_points = object@d_obs,
    ...
  )

  if (S7::S7_inherits(vario, UtopVariogramCloud)) {
    object@variogram_cloud <- vario
  } else {
    object@variogram <- vario
  }

  object
}

S7::method(utop_variogram, utop_sf_class) <- function(
  object,
  formula = NULL,
  params = NULL,
  cloud = NULL,
  disc_points = NULL,
  ...
) {
  if (is.null(params)) {
    params <- UtopParams()
  }
  if (!S7::S7_inherits(params, UtopParams)) {
    params <- utop_params(params, ...)
  }
  if (is.null(cloud)) {
    cloud <- params@cloud
  }
  compute_variogram_sf(
    object = object,
    formula = formula,
    params = params,
    cloud = cloud,
    disc_points = disc_points,
    ...
  )
}

S7::method(utop_variogram, utop_stars_class) <- function(
  object,
  formula = NULL,
  params = NULL,
  cloud = NULL,
  abins = NULL,
  dbins = NULL,
  disc_points = NULL,
  ...
) {
  if (is.null(params)) {
    params <- UtopParams()
  }
  if (!S7::S7_inherits(params, UtopParams)) {
    params <- utop_params(params, ...)
  }
  if (is.null(cloud)) {
    cloud <- params@cloud
  }
  compute_variogram_stars(
    object = object,
    formula = formula,
    params = params,
    cloud = cloud,
    abins = abins,
    dbins = dbins,
    disc_points = disc_points,
    ...
  )
}
