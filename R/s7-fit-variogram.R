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
fit_variogram_impl <- function(state, params = NULL, iprint = 0, ...) {
  state$params <- utop_replace_params(
    current = state$params,
    params = params,
    observations = state$observations,
    formula = state$formula
  )
  params_obj <- state$params

  if (
    (params_obj@cloud && is.null(state$variogram_cloud)) ||
      (!params_obj@cloud && is.null(state$variogram))
  ) {
    obj <- utop_variogram(utop_utop_from_state(state), ...)
    state <- utop_state_from_utop(obj)
  }

  if (params_obj@nugget && is.null(state$overlap_obs)) {
    state$overlap_obs <- find_overlap(
      state$observations,
      state$observations,
      partial_overlap = TRUE
    )
  }

  if (params_obj@cloud) {
    dists <- if (params_obj@g_dist_est && !is.null(state$g_dist_obs)) {
      state$g_dist_obs
    } else if (!is.null(state$d_obs)) {
      state$d_obs
    } else {
      NULL
    }
    a_over <- if (params_obj@nugget) state$overlap_obs else NULL
    fit <- fit_variogram_cloud(
      state$variogram_cloud,
      observations = state$observations,
      params = params_obj,
      dists = dists,
      a_over = a_over,
      iprint = iprint,
      ...
    )
  } else {
    dists <- if (params_obj@g_dist_est && !is.null(state$g_dist_bin)) {
      state$g_dist_bin
    } else if (!is.null(state$d_bin)) {
      state$d_bin
    } else {
      NULL
    }
    a_over <- if (params_obj@nugget) {
      findVarioOverlap(utop_variogram_data(state$variogram))
    } else {
      NULL
    }
    fit <- fit_variogram_binned(
      state$variogram,
      params = params_obj,
      dists = dists,
      a_over = a_over,
      iprint = iprint,
      ...
    )
  }

  state$variogram_model <- fit$model
  state$var_fit <- fit$var_fit
  if (!is.null(fit$caches$d_bin)) {
    state$d_bin <- fit$caches$d_bin
  }
  if (!is.null(fit$caches$g_dist_bin)) {
    state$g_dist_bin <- fit$caches$g_dist_bin
  }
  if (!is.null(fit$caches$d_obs)) {
    state$d_obs <- fit$caches$d_obs
  }
  if (!is.null(fit$caches$g_dist_obs)) {
    state$g_dist_obs <- fit$caches$g_dist_obs
  }

  state
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
  state <- fit_variogram_impl(
    utop_state_from_utop(object),
    params = params,
    iprint = iprint,
    ...
  )
  utop_utop_from_state(state)
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
