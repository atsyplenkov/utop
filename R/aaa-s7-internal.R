#' @noRd
utop_param_names <- c(
  model = "model",
  nugget = "nugget",
  unc = "unc",
  r_resol = "rresol",
  h_resol = "hresol",
  rs_type = "rstype",
  hs_type = "hstype",
  cloud = "cloud",
  amul = "amul",
  dmul = "dmul",
  fit_method = "fit.method",
  g_dist_est = "gDistEst",
  g_dist_pred = "gDistPred",
  var_clean = "varClean",
  max_dist = "maxdist",
  n_max = "nmax",
  n_clus = "nclus",
  cn_areas = "cnAreas",
  clus_type = "clusType",
  outfile = "outfile",
  partial_overlap = "partialOverlap",
  wlim = "wlim",
  wlim_method = "wlimMethod",
  singular_solve = "singularSolve",
  uk_trend_support = "ukTrendSupport",
  cv = "cv",
  debug_level = "debug.level",
  par_init = "parInit"
)

#' @noRd
coerce_utop_params <- function(params = NULL, ...) {
  dots <- list(...)
  if ("newPar" %in% names(dots)) {
    new_params <- dots$newPar
    dots$newPar <- NULL
    base <- if (S7::S7_inherits(params, UtopParams)) {
      params
    } else {
      utop_params(params = params, ...)
    }
    return(utop_params(params = base, new_params = new_params, ...))
  }

  if (S7::S7_inherits(params, UtopParams)) {
    if (length(dots) > 0L) {
      param_names <- c(
        names(S7::props(params)),
        unname(utop_param_names),
        "gDist"
      )
      param_dots <- dots[intersect(names(dots), param_names)]
      if (length(param_dots) > 0L) {
        return(utop_params(
          params = params,
          new_params = utop_params_apply_aliases(param_dots)
        ))
      }
    }
    return(params)
  }

  if (is.list(params) && length(params) > 0L) {
    param_names <- c(
      names(S7::props(UtopParams())),
      unname(utop_param_names),
      "gDist"
    )
    param_dots <- dots[intersect(names(dots), param_names)]
    if (length(param_dots) > 0L) {
      return(utop_params(
        params = params,
        new_params = utop_params_apply_aliases(param_dots)
      ))
    }
    return(utop_params(params = params))
  }

  if (is.null(params)) {
    return(utop_params(...))
  }

  utop_params(params = params, ...)
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
utop_variogram_model_from_rtop <- function(model) {
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

  unknown <- setdiff(names(values), names(S7::props(params)))
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
utop_params_apply_aliases <- function(values) {
  if (is.null(values) || length(values) == 0L) {
    return(values)
  }
  if ("geoDist" %in% names(values)) {
    stop("geoDist is not used anymore, please use g_dist", call. = FALSE)
  }
  if ("gDist" %in% names(values)) {
    g_dist <- isTRUE(values$gDist)
    values$g_dist_est <- g_dist
    values$g_dist_pred <- g_dist
    values$gDist <- NULL
  }

  legacy_map <- stats::setNames(
    names(utop_param_names),
    unname(utop_param_names)
  )
  legacy_names <- intersect(names(values), names(legacy_map))
  for (old_name in legacy_names) {
    new_name <- legacy_map[[old_name]]
    if (old_name == new_name) {
      next
    }
    if (is.null(values[[new_name]])) {
      values[[new_name]] <- values[[old_name]]
    }
    values[[old_name]] <- NULL
  }

  values
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
