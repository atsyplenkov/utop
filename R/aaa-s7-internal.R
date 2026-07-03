#' @noRd
utop_param_fields <- function() {
  names(S7::props(UtopParams()))
}

#' @noRd
utop_require_named_values <- function(values, caller = "utop_params()") {
  if (length(values) == 0L) {
    return(values)
  }
  if (
    is.null(names(values)) || anyNA(names(values)) || any(names(values) == "")
  ) {
    stop(caller, " only accepts named arguments", call. = FALSE)
  }
  values
}

#' @noRd
utop_is_default_par_init <- function(par_init, model) {
  if (is.null(par_init)) {
    return(TRUE)
  }
  default <- find_par_init_default(model)
  isTRUE(all.equal(par_init, default, check.attributes = FALSE))
}

#' @noRd
utop_complete_params <- function(params, observations = NULL, formula = NULL) {
  if (!S7::S7_inherits(params, UtopParams)) {
    stop("internal error: expected a UtopParams object", call. = FALSE)
  }

  if (
    !is.null(observations) &&
      utop_is_default_par_init(params@par_init, params@model)
  ) {
    if (is.null(formula)) {
      formula <- as.formula(utop_default_formula(observations))
    }
    params@par_init <- find_par_init(
      formula = formula,
      observations = observations,
      model = params@model
    )
  } else if (is.null(params@par_init)) {
    params@par_init <- find_par_init_default(params@model)
  }

  params
}

#' @noRd
utop_require_params <- function(
  params = NULL,
  observations = NULL,
  formula = NULL,
  arg = "params"
) {
  if (is.null(params)) {
    return(utop_complete_params(
      UtopParams(),
      observations = observations,
      formula = formula
    ))
  }
  if (!S7::S7_inherits(params, UtopParams)) {
    stop(
      "`",
      arg,
      "` must be a UtopParams object. Create it with utop_params().",
      call. = FALSE
    )
  }

  utop_complete_params(params, observations = observations, formula = formula)
}

#' @noRd
utop_replace_params <- function(
  current = NULL,
  params = NULL,
  observations = NULL,
  formula = NULL,
  arg = "params"
) {
  if (is.null(params)) {
    if (is.null(current)) {
      return(utop_require_params(
        NULL,
        observations = observations,
        formula = formula,
        arg = arg
      ))
    }
    return(utop_complete_params(
      current,
      observations = observations,
      formula = formula
    ))
  }

  utop_require_params(
    params,
    observations = observations,
    formula = formula,
    arg = arg
  )
}

#' @noRd
coerce_variogram_model <- function(model) {
  if (is.null(model)) {
    return(NULL)
  }
  if (S7::S7_inherits(model, UtopVariogramModel)) {
    out <- list(model = model@model, params = model@params)
    class(out) <- "rtopVariogramModel"
    if (!is.null(model@ss_err)) {
      attr(out, "SSErr") <- model@ss_err
    }
    if (!is.null(model@criterion)) {
      attr(out, "criterion") <- model@criterion
    }
    return(out)
  }
  model
}

#' @noRd
variogram_is_cloud <- function(x) {
  if (S7::S7_inherits(x, UtopVariogramCloud)) {
    return(TRUE)
  }
  if (S7::S7_inherits(x, UtopVariogram)) {
    return(FALSE)
  }
  inherits(x, "data.frame") && "ord" %in% names(x)
}

#' @noRd
utop_tag_variogram_class <- function(data, cloud = FALSE) {
  class(data) <- c(
    if (cloud) "rtopVariogramCloud" else "rtopVariogram",
    "data.frame"
  )
  data
}

#' @noRd
utop_variogram_model_from_list <- function(model) {
  if (is.null(model)) {
    return(NULL)
  }
  if (S7::S7_inherits(model, UtopVariogramModel)) {
    return(model)
  }
  UtopVariogramModel(
    model = model$model,
    params = model$params,
    ss_err = attr(model, "SSErr", exact = TRUE),
    criterion = attr(model, "criterion", exact = TRUE)
  )
}

#' @noRd
update_utop_params <- function(params, values) {
  if (is.null(values) || length(values) == 0L) {
    return(params)
  }

  values <- utop_require_named_values(values)
  unknown <- setdiff(names(values), utop_param_fields())
  if (length(unknown) > 0L) {
    stop(
      "unknown utop parameter(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  for (name in names(values)) {
    value <- values[[name]]
    if (!is.null(value)) {
      S7::prop(params, name) <- value
    }
  }

  params
}

#' @noRd
utop_wrap_variogram_result <- function(data, cloud = FALSE) {
  if (cloud) {
    return(UtopVariogramCloud(data = data))
  }
  UtopVariogram(data = data)
}

#' @noRd
utop_variogram_data <- function(object) {
  if (S7::S7_inherits(object, UtopVariogram)) {
    return(object@data)
  }
  if (S7::S7_inherits(object, UtopVariogramCloud)) {
    return(object@data)
  }
  stop("expected a UtopVariogram or UtopVariogramCloud object", call. = FALSE)
}

#' @noRd
utop_support_bbox <- function(observations, prediction_locations = NULL) {
  if (is.null(prediction_locations)) {
    return(utop_bbox(observations))
  }

  obs_bbox <- utop_bbox(observations)
  pred_bbox <- utop_bbox(prediction_locations)

  c(
    xmin = min(obs_bbox[[1]], pred_bbox[[1]]),
    ymin = min(obs_bbox[[2]], pred_bbox[[2]]),
    xmax = max(obs_bbox[[3]], pred_bbox[[3]]),
    ymax = max(obs_bbox[[4]], pred_bbox[[4]])
  )
}

#' @noRd
check_matching_crs <- function(observations, prediction_locations) {
  obs_crs <- sf::st_crs(utop_as_sf(observations))
  pred_crs <- sf::st_crs(utop_as_sf(prediction_locations))
  if (!isTRUE(all.equal(is.na(obs_crs), is.na(pred_crs)))) {
    stop(
      "only one of observations and prediction_locations have projection",
      call. = FALSE
    )
  }
  if (!is.na(obs_crs) && obs_crs != pred_crs) {
    stop(
      "observations and prediction_locations have different projections: ",
      obs_crs,
      " ",
      pred_crs,
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' @noRd
utop_state_has <- function(state, name) {
  !is.null(state[[name]])
}

#' @noRd
utop_utop_from_state <- function(state) {
  args <- list(params = state$params)
  slots <- c(
    "observations",
    "prediction_locations",
    "formula",
    "variogram",
    "variogram_cloud",
    "variogram_model",
    "predictions",
    "simulations",
    "d_obs",
    "d_pred",
    "g_dist_obs",
    "g_dist_pred",
    "g_dist_pred_obs",
    "g_dist_bin",
    "d_bin",
    "overlap_obs",
    "overlap_pred_obs",
    "var_mat_obs",
    "var_mat_pred_obs",
    "var_mat_pred",
    "var_fit",
    "cv_info",
    "weight",
    "removed",
    "uk_residual",
    "check_vario"
  )
  for (name in slots) {
    if (utop_state_has(state, name)) {
      args[[name]] <- state[[name]]
    }
  }
  do.call(Utop, args)
}

#' @noRd
utop_state_from_utop <- function(object) {
  list(
    observations = object@observations,
    prediction_locations = object@prediction_locations,
    formula = object@formula,
    params = object@params,
    variogram = object@variogram,
    variogram_cloud = object@variogram_cloud,
    variogram_model = object@variogram_model,
    predictions = object@predictions,
    simulations = object@simulations,
    d_obs = object@d_obs,
    d_pred = object@d_pred,
    g_dist_obs = object@g_dist_obs,
    g_dist_pred = object@g_dist_pred,
    g_dist_pred_obs = object@g_dist_pred_obs,
    g_dist_bin = object@g_dist_bin,
    d_bin = object@d_bin,
    overlap_obs = object@overlap_obs,
    overlap_pred_obs = object@overlap_pred_obs,
    var_mat_obs = object@var_mat_obs,
    var_mat_pred_obs = object@var_mat_pred_obs,
    var_mat_pred = object@var_mat_pred,
    var_fit = object@var_fit,
    cv_info = object@cv_info,
    weight = object@weight,
    removed = object@removed,
    uk_residual = object@uk_residual,
    check_vario = object@check_vario
  )
}
