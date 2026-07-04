#' @noRd
fit_variogram_binned <- function(
  variogram,
  params,
  dists = NULL,
  a_over = NULL,
  iprint = 0,
  ...
) {
  vario <- utop_variogram_data(variogram)
  caches <- list()

  if (is.null(dists)) {
    dists <- utop_disc(variogram, params = params)
    caches$d_bin <- dists
  }
  if (params@g_dist_est && !is.matrix(dists)) {
    dists <- compute_g_dist_list(dists, params = params, ...)
    caches$g_dist_bin <- dists
  }
  if (params@nugget && is.null(a_over)) {
    a_over <- findVarioOverlap(vario)
  }

  implicit <- utop_model_implicit(params@model)

  scres <- sceua::sceua(
    objfunc,
    params@par_init[, 3],
    lower = params@par_init[, 1],
    upper = params@par_init[, 2],
    varioIn = vario,
    dists = dists,
    aOver = a_over,
    gDistEst = params@g_dist_est,
    model = params@model,
    resol = params@h_resol,
    fit.method = params@fit_method,
    implicit = implicit,
    iprint = iprint,
    cloud = FALSE,
    ...
  )
  best_par <- scres$par
  fit <- scres$value
  vf <- objfunc(
    best_par,
    varioIn = vario,
    dists = dists,
    aOver = a_over,
    gDistEst = params@g_dist_est,
    last = TRUE,
    model = params@model,
    resol = params@h_resol,
    cloud = FALSE,
    ...
  )

  model <- UtopVariogramModel(
    model = params@model,
    params = best_par,
    ss_err = vf$errSum,
    criterion = fit
  )
  if (!params@nugget) {
    model@params[3] <- 0
  }

  list(model = model, var_fit = vf$varFit, caches = caches)
}

#' @noRd
fit_variogram_cloud <- function(
  variogram,
  observations,
  params,
  dists = NULL,
  a_over = NULL,
  iprint = 0,
  ...
) {
  vario <- utop_variogram_data(variogram)
  caches <- list()

  if (is.null(dists)) {
    dists <- utop_disc(observations, params = params, ...)
    caches$d_obs <- dists
  }
  if (params@g_dist_est && is.list(dists)) {
    dists <- compute_g_dist_list(dists, params = params, ...)
    caches$g_dist_obs <- dists
  }

  implicit <- utop_model_implicit(params@model)

  scres <- sceua::sceua(
    objfunc,
    params@par_init[, 3],
    lower = params@par_init[, 1],
    upper = params@par_init[, 2],
    varioIn = vario,
    dists = dists,
    aOver = a_over,
    gDist = params@g_dist_est,
    model = params@model,
    fit.method = params@fit_method,
    debug.level = params@debug_level,
    iprint = iprint,
    cloud = TRUE,
    ...
  )
  best_par <- scres$par
  vf <- objfunc(
    best_par,
    varioIn = vario,
    dists = dists,
    aOver = a_over,
    gDistEst = params@g_dist_est,
    last = TRUE,
    model = params@model,
    debug.level = params@debug_level,
    cloud = TRUE,
    ...
  )

  model <- UtopVariogramModel(
    model = params@model,
    params = best_par,
    ss_err = vf$errSum
  )
  if (!params@nugget) {
    model@params[3] <- 0
  }

  list(model = model, var_fit = vf$varFit, caches = caches)
}

#' @noRd
fit_variogram_impl <- function(object, params = NULL, iprint = 0, ...) {
  object@params <- utop_replace_params(
    current = object@params,
    params = params,
    observations = object@observations,
    formula = object@formula
  )
  params_obj <- object@params

  if (
    (params_obj@cloud && is.null(object@variogram_cloud)) ||
      (!params_obj@cloud && is.null(object@variogram))
  ) {
    object <- utop_variogram(object, ...)
  }

  if (params_obj@nugget && is.null(object@overlap_obs)) {
    object@overlap_obs <- find_overlap(
      object@observations,
      object@observations,
      partial_overlap = TRUE
    )
  }

  if (params_obj@cloud) {
    dists <- if (params_obj@g_dist_est && !is.null(object@g_dist_obs)) {
      object@g_dist_obs
    } else if (!is.null(object@d_obs)) {
      object@d_obs
    } else {
      NULL
    }
    a_over <- if (params_obj@nugget) object@overlap_obs else NULL
    fit <- fit_variogram_cloud(
      object@variogram_cloud,
      observations = object@observations,
      params = params_obj,
      dists = dists,
      a_over = a_over,
      iprint = iprint,
      ...
    )
  } else {
    dists <- if (params_obj@g_dist_est && !is.null(object@g_dist_bin)) {
      object@g_dist_bin
    } else if (!is.null(object@d_bin)) {
      object@d_bin
    } else {
      NULL
    }
    a_over <- if (params_obj@nugget) {
      findVarioOverlap(utop_variogram_data(object@variogram))
    } else {
      NULL
    }
    fit <- fit_variogram_binned(
      object@variogram,
      params = params_obj,
      dists = dists,
      a_over = a_over,
      iprint = iprint,
      ...
    )
  }

  object@variogram_model <- fit$model
  object@var_fit <- fit$var_fit
  if (!is.null(fit$caches$d_bin)) {
    object@d_bin <- fit$caches$d_bin
  }
  if (!is.null(fit$caches$g_dist_bin)) {
    object@g_dist_bin <- fit$caches$g_dist_bin
  }
  if (!is.null(fit$caches$d_obs)) {
    object@d_obs <- fit$caches$d_obs
  }
  if (!is.null(fit$caches$g_dist_obs)) {
    object@g_dist_obs <- fit$caches$g_dist_obs
  }

  object
}

#' Fit a variogram model
#'
#' @param object A [Utop], [UtopVariogram], or [UtopVariogramCloud] object.
#' @param ... Method-specific arguments. For [Utop] objects, `iprint`
#'   controls optimiser verbosity.
#'
#' @return A fitted object.
#' @export
utop_fit_variogram <- S7::new_generic(
  name = "utop_fit_variogram",
  dispatch_args = "object",
  fun = function(object, ...) {
    S7::S7_dispatch()
  }
)

S7::method(utop_fit_variogram, Utop) <- function(
  object,
  params = NULL,
  iprint = 0,
  ...
) {
  fit_variogram_impl(object, params = params, iprint = iprint, ...)
}

S7::method(utop_fit_variogram, UtopVariogram) <- function(
  object,
  params = NULL,
  iprint = 0,
  ...
) {
  params <- utop_require_params(params)
  fit_variogram_binned(object, params = params, iprint = iprint, ...)$model
}

S7::method(utop_fit_variogram, UtopVariogramCloud) <- function(
  object,
  observations,
  params = NULL,
  iprint = 0,
  ...
) {
  params <- utop_require_params(params, observations = observations)
  fit_variogram_cloud(
    object,
    observations = observations,
    params = params,
    iprint = iprint,
    ...
  )$model
}

S7::method(utop_fit_variogram, utop_sf_class) <- function(
  object,
  params = NULL,
  iprint = 0,
  ...
) {
  params <- utop_require_params(params, observations = object)
  vario <- utop_variogram(object, params = params, ...)
  utop_fit_variogram(vario, params = params, iprint = iprint, ...)
}

S7::method(utop_fit_variogram, utop_stars_class) <- function(
  object,
  params = NULL,
  iprint = 0,
  ...
) {
  params <- utop_require_params(params, observations = object)
  vario <- utop_variogram(object, params = params, ...)
  utop_fit_variogram(vario, params = params, iprint = iprint, ...)
}
