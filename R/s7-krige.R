#' @noRd
compute_krige_utop <- function(
  object,
  var_mat_update = FALSE,
  params = list(),
  ...
) {
  split <- utop_krige_extra_params(params, ...)
  object@params <- utop_params(object@params, new_params = split$params)
  params_for_compute <- utop_krige_apply_runtime(object@params, split$runtime)
  runtime <- utop_krige_runtime_only(split$runtime)

  if (!is.null(runtime$nsim) && runtime$nsim > 0) {
    return(compute_sim_utop(
      object,
      var_mat_update = var_mat_update,
      params = params,
      ...
    ))
  }

  if (
    is.null(object@var_mat_obs) ||
      is.null(object@var_mat_pred_obs) ||
      var_mat_update
  ) {
    object <- utop_var_mat(
      object,
      var_mat_update = var_mat_update,
      params = split$params,
      split$dots
    )
  }

  krige_args <- c(
    list(
      object = object@observations,
      predictionLocations = object@prediction_locations,
      varMatObs = object@var_mat_obs,
      varMatPredObs = object@var_mat_pred_obs,
      params = params_for_compute,
      formulaString = object@formula,
      dObs = object@d_obs,
      dPred = object@d_pred
    ),
    runtime,
    split$dots
  )

  krige_res <- if (inherits(object@observations, "stars")) {
    do.call(compute_krige_stars, krige_args)
  } else {
    do.call(compute_krige, krige_args)
  }

  object@predictions <- krige_res$predictions
  if ("cvInfo" %in% names(krige_res)) {
    object@cv_info <- krige_res$cvInfo
  }
  if ("weight" %in% names(krige_res)) {
    object@weight <- krige_res$weight
  }
  if ("removed" %in% names(krige_res)) {
    object@removed <- krige_res$removed
  }
  object
}

#' Krige with a utop object
#'
#' @param object A fitted [Utop] object or spatial input.
#' @param ... Method-specific arguments. For [Utop] objects, `var_mat_update`
#'   and `params` control matrix recomputation and parameter updates.
#'
#' @return A [Utop] object with predictions attached, or an `sf` result.
#' @export
utop_krige <- S7::new_generic(
  name = "utop_krige",
  dispatch_args = "object",
  fun = function(object, ...) {
    S7::S7_dispatch()
  }
)

S7::method(utop_krige, Utop) <- function(
  object,
  var_mat_update = FALSE,
  params = list(),
  ...
) {
  compute_krige_utop(
    object,
    var_mat_update = var_mat_update,
    params = params,
    ...
  )
}

#' @noRd
utop_krige_apply_runtime <- function(params, runtime) {
  if (length(runtime) == 0L) {
    return(params)
  }
  runtime <- utop_params_apply_aliases(runtime)
  updates <- runtime[intersect(names(runtime), names(S7::props(params)))]
  if (length(updates) > 0L) {
    params <- update_utop_params(params, updates)
  }
  params
}

#' @noRd
utop_krige_runtime_only <- function(runtime) {
  krige_only <- c(
    "nsim",
    "lambda",
    "olags",
    "plags",
    "lagExact",
    "sel",
    "wret"
  )
  runtime[intersect(names(runtime), krige_only)]
}

#' @noRd
utop_krige_extra_params <- function(params, ...) {
  dots <- list(...)
  runtime_names <- c(
    "cv",
    "wlim",
    "wlimMethod",
    "wlim_method",
    "nmax",
    "n_max",
    "maxdist",
    "max_dist",
    "nsim",
    "lambda",
    "singularSolve",
    "singular_solve",
    "olags",
    "plags",
    "lagExact",
    "sel",
    "wret"
  )
  if (is.null(params)) {
    params <- list()
  }
  runtime <- c(params, dots)
  runtime <- runtime[intersect(names(runtime), runtime_names)]
  utop_updates <- params[!names(params) %in% runtime_names]
  dots <- dots[!names(dots) %in% runtime_names]
  list(
    params = utop_updates,
    runtime = runtime,
    dots = dots
  )
}

S7::method(utop_krige, utop_sf_class) <- function(
  object,
  prediction_locations = NULL,
  params = NULL,
  formula = NULL,
  var_mat_obs = NULL,
  var_mat_pred_obs = NULL,
  ...
) {
  split <- utop_krige_extra_params(params, ...)
  params_obj <- if (is.null(split$params)) {
    UtopParams()
  } else if (S7::S7_inherits(split$params, UtopParams)) {
    split$params
  } else {
    utop_params(split$params)
  }
  utop_obj <- utop_object(
    observations = object,
    prediction_locations = prediction_locations,
    formula = formula,
    params = params_obj
  )
  if (!is.null(var_mat_obs)) {
    utop_obj@var_mat_obs <- var_mat_obs
    model <- attr(var_mat_obs, "variogramModel", exact = TRUE)
    if (!is.null(model) && is.null(utop_obj@variogram_model)) {
      utop_obj@variogram_model <- utop_variogram_model_from_rtop(model)
    }
  }
  if (!is.null(var_mat_pred_obs)) {
    utop_obj@var_mat_pred_obs <- var_mat_pred_obs
  }
  do.call(utop_krige, c(list(utop_obj), split$runtime, split$dots))
}

S7::method(utop_krige, utop_stars_class) <- function(
  object,
  prediction_locations = NULL,
  params = NULL,
  formula = NULL,
  var_mat_obs = NULL,
  var_mat_pred_obs = NULL,
  ...
) {
  split <- utop_krige_extra_params(params, ...)
  params_obj <- if (is.null(split$params)) {
    UtopParams()
  } else if (S7::S7_inherits(split$params, UtopParams)) {
    split$params
  } else {
    utop_params(split$params)
  }
  utop_obj <- utop_object(
    observations = object,
    prediction_locations = prediction_locations,
    formula = formula,
    params = params_obj
  )
  if (!is.null(var_mat_obs)) {
    utop_obj@var_mat_obs <- var_mat_obs
    model <- attr(var_mat_obs, "variogramModel", exact = TRUE)
    if (!is.null(model) && is.null(utop_obj@variogram_model)) {
      utop_obj@variogram_model <- utop_variogram_model_from_rtop(model)
    }
  }
  if (!is.null(var_mat_pred_obs)) {
    utop_obj@var_mat_pred_obs <- var_mat_pred_obs
  }
  do.call(utop_krige, c(list(utop_obj), split$runtime, split$dots))
}
