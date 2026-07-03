#' @noRd
compute_krige_utop <- function(
  object,
  var_mat_update = FALSE,
  params = NULL,
  ...
) {
  split <- utop_krige_extra_params(params, ...)
  object@params <- utop_replace_params(
    current = object@params,
    params = split$params,
    observations = object@observations,
    formula = object@formula
  )
  params_for_compute <- utop_krige_apply_runtime(object@params, split$runtime)
  runtime <- utop_krige_runtime_only(split$runtime)

  if (!is.null(runtime$nsim) && runtime$nsim > 0) {
    return(do.call(
      compute_sim_utop,
      c(
        list(
          object = object,
          var_mat_update = var_mat_update,
          params = split$params
        ),
        split$runtime,
        split$dots
      )
    ))
  }

  if (
    is.null(object@var_mat_obs) ||
      is.null(object@var_mat_pred_obs) ||
      var_mat_update
  ) {
    object <- do.call(
      utop_var_mat,
      c(
        list(
          object = object,
          var_mat_update = var_mat_update,
          params = object@params
        ),
        split$dots
      )
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
  params = NULL,
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
  updates <- runtime[intersect(names(runtime), utop_param_fields())]
  if (length(updates) > 0L) {
    params <- update_utop_params(params, updates)
  }
  params
}

#' @noRd
utop_krige_runtime_only <- function(runtime) {
  if ("lag_exact" %in% names(runtime)) {
    runtime$lagExact <- runtime$lag_exact
    runtime$lag_exact <- NULL
  }
  krige_only <- c("nsim", "lambda", "olags", "plags", "lagExact", "sel", "wret")
  runtime[intersect(names(runtime), krige_only)]
}

#' @noRd
utop_krige_extra_params <- function(params, ...) {
  dots <- list(...)
  runtime_names <- c(
    "cv",
    "nsim",
    "lambda",
    "olags",
    "plags",
    "lag_exact",
    "sel",
    "wret"
  )
  runtime <- dots[intersect(names(dots), runtime_names)]
  dots <- dots[!names(dots) %in% runtime_names]

  inline_params <- intersect(names(dots), utop_param_fields())
  if (length(inline_params) > 0L) {
    stop(
      "pass utop parameters via a UtopParams object in `params`, not as individual arguments: ",
      paste(inline_params, collapse = ", "),
      call. = FALSE
    )
  }

  list(params = params, runtime = runtime, dots = dots)
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
  params_obj <- utop_require_params(
    split$params,
    observations = object,
    formula = formula
  )
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
      utop_obj@variogram_model <- utop_variogram_model_from_list(model)
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
  params_obj <- utop_require_params(
    split$params,
    observations = object,
    formula = formula
  )
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
      utop_obj@variogram_model <- utop_variogram_model_from_list(model)
    }
  }
  if (!is.null(var_mat_pred_obs)) {
    utop_obj@var_mat_pred_obs <- var_mat_pred_obs
  }
  do.call(utop_krige, c(list(utop_obj), split$runtime, split$dots))
}
