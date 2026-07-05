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

# Split `...` into (params, runtime, dots) for functions that accept a mix of
# runtime-only arguments (e.g. nsim, cv) and passthrough dots. Any utop
# parameter name found outside `params` is rejected with an error pointing the
# caller at the UtopParams API. `runtime_names` is the per-call site list of
# argument names that should be peeled off as runtime arguments.
#' @noRd
utop_split_runtime <- function(params, runtime_names, ...) {
  dots <- list(...)
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

# Coerce any variogram-like input (S7 UtopVariogram[UtopVariogramCloud], a
# tagged data.frame, or a plain data.frame) to a plain data.frame. Used at
# boundaries that consume variograms with S3 `$` / `[` access, which S7 objects
# do not support. NULL passes through unchanged.
#' @noRd
utop_as_variogram_df <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (S7::S7_inherits(x, UtopVariogram) || S7::S7_inherits(x, UtopVariogramCloud)) {
    x <- x@data
  }
  as.data.frame(x)
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
