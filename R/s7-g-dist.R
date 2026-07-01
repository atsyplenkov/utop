#' @noRd
compute_g_dist_list <- function(
  object,
  object2 = NULL,
  diag = FALSE,
  debug.level = 0,
  params = NULL,
  ...
) {
  variogram_model <- list(model = "Gho", params = 0)
  if (debug.level == 1) {
    print("Creating Ghos distances. This can take some time")
  }
  if (inherits(object[[1]], "sf") || inherits(object[[1]], "sfc")) {
    g_dist <- compute_var_mat_list(
      object,
      object2,
      diag = diag,
      variogramModel = variogram_model,
      debug.level = debug.level,
      params = params,
      ...
    )
  } else {
    g_dist <- data.frame(c1 = c(rep(0, length(object))), c2 = 0, cb = 0)
    l_areas <- lapply(object, FUN = function(aa) aa[[1]])
    g_dist[, 1] <- mapply(
      vred,
      l_areas,
      MoreArgs = list(vredTyp = "ind", variogramModel = variogram_model)
    )
    l_areas <- lapply(object, FUN = function(aa) aa[[2]])
    g_dist[, 2] <- mapply(
      vred,
      l_areas,
      MoreArgs = list(vredTyp = "ind", variogramModel = variogram_model)
    )
    g_dist[, 3] <- mapply(
      vred,
      object,
      MoreArgs = list(vredTyp = "ind", variogramModel = variogram_model)
    )
  }
  as.matrix(g_dist)
}

#' @noRd
compute_g_dist_sf <- function(object, object2 = NULL, ...) {
  params <- UtopParams()
  d_obs <- utop_disc(object, params = params, ...)
  g_dist_obs <- compute_g_dist_list(d_obs, params = params, ...)
  if (!is.null(object2)) {
    d_pred <- utop_disc(object2, params = params, ...)
    g_dist_pred_obs <- compute_g_dist_list(d_obs, d_pred, params = params, ...)
    g_dist_pred <- compute_g_dist_list(
      d_pred,
      diag = TRUE,
      params = params,
      ...
    )
    return(list(
      g_dist_obs = g_dist_obs,
      g_dist_pred = g_dist_pred,
      g_dist_pred_obs = g_dist_pred_obs
    ))
  }
  list(g_dist_obs = g_dist_obs)
}

#' Compute geostatistical distances
#'
#' @param object A [Utop] object, spatial object, or discretisation list.
#' @param ... Method-specific arguments. For [Utop] objects, `params` updates
#'   parameters.
#'
#' @return Geostatistical distances.
#' @export
utop_g_dist <- S7::new_generic(
  name = "utop_g_dist",
  dispatch_args = "object",
  fun = function(object, ...) {
    S7::S7_dispatch()
  }
)

S7::method(utop_g_dist, Utop) <- function(object, params = list(), ...) {
  object@params <- utop_params(object@params, new_params = params, ...)
  debug.level <- if (object@params@debug_level > 1) {
    object@params@debug_level
  } else {
    0
  }

  if (is.null(object@d_obs)) {
    object@d_obs <- utop_disc(object@observations, params = object@params, ...)
  }
  if (
    is.null(object@d_pred) &&
      !is.null(object@prediction_locations)
  ) {
    object@d_pred <- utop_disc(
      object@prediction_locations,
      params = object@params,
      ...
    )
  }

  object@g_dist_obs <- compute_g_dist_list(
    object@d_obs,
    object@d_obs,
    debug.level = debug.level,
    params = object@params,
    ...
  )

  if (!is.null(object@d_pred)) {
    object@g_dist_pred_obs <- compute_g_dist_list(
      object@d_obs,
      object@d_pred,
      debug.level = debug.level,
      params = object@params,
      ...
    )
    object@g_dist_pred <- compute_g_dist_list(
      object@d_pred,
      object@d_pred,
      diag = TRUE,
      debug.level = debug.level,
      params = object@params,
      ...
    )
  }

  object
}

S7::method(utop_g_dist, utop_sf_class) <- function(
  object,
  object2 = NULL,
  ...
) {
  compute_g_dist_sf(object, object2 = object2, ...)
}

S7::method(utop_g_dist, utop_stars_class) <- function(
  object,
  object2 = NULL,
  ...
) {
  object <- utop_stars_support(object)
  if (!is.null(object2)) {
    object2 <- utop_as_sf(object2)
  }
  compute_g_dist_sf(object, object2 = object2, ...)
}

S7::method(utop_g_dist, S7::class_list) <- function(object, object2 = NULL, ...) {
  compute_g_dist_list(object, object2 = object2, ...)
}
