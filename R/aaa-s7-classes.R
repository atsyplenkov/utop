#' S7 classes for utop objects
#'
#' These classes provide the formal object model used by the snake_case utop
#' API. Prefer `utop_object()` and related helpers for user-facing workflows.
#'
#' @param model Variogram model name.
#' @param nugget,unc,cloud,cv Logical parameters.
#' @param r_resol,h_resol,amul,dmul,n_max,n_clus Numeric or integer tuning
#'   parameters.
#' @param cn_areas Number of areas used in cluster-related calculations.
#' @param rs_type,hs_type,clus_type,outfile,wlim_method,uk_trend_support
#'   Character parameters.
#' @param fit_method,g_dist_est,g_dist_pred,var_clean,max_dist Variogram and
#'   distance parameters.
#' @param partial_overlap Whether partially overlapping areas are supported.
#' @param wlim,singular_solve,debug_level Additional tuning parameters.
#' @param par_init Initial parameter bounds.
#' @param data Data stored in variogram objects.
#' @param params Numeric variogram model parameters.
#' @param ss_err,criterion Variogram fitting diagnostics.
#' @param observations,prediction_locations Spatial observations and prediction
#'   locations.
#' @param formula Model formula.
#' @param variogram,variogram_cloud,variogram_model Variogram results.
#' @param predictions,simulations Prediction and simulation results.
#' @param d_obs,d_pred,g_dist_obs,g_dist_pred,g_dist_pred_obs,g_dist_bin,d_bin
#'   Discretisation and geostatistical distance objects.
#' @param overlap_obs,overlap_pred_obs Overlap matrices.
#' @param var_mat_obs,var_mat_pred_obs,var_mat_pred Covariance matrices.
#' @param var_fit,cv_info,weight,removed,uk_residual,check_vario Diagnostics
#'   and intermediate results.
#' @param method_parameters Method parameter summary.
#'
#' @name utop-s7-classes
#' @aliases Utop
#' @aliases UtopParams
#' @aliases UtopVariogram
#' @aliases UtopVariogramCloud
#' @aliases UtopVariogramModel
#' @rawNamespace if (getRversion() < "4.3.0") importFrom("S7", "@")
#' @importFrom S7 S7_dispatch S7_inherits class_any class_character
#' @importFrom S7 class_data.frame class_formula class_logical class_numeric
#' @importFrom S7 method methods_register new_class new_generic
#' @importFrom S7 new_property prop
NULL

utop_optional <- function(class = S7::class_any) {
  S7::new_property(NULL | class, default = NULL)
}

#' @rdname utop-s7-classes
#' @export
UtopParams <- S7::new_class(
  "UtopParams",
  package = "utop",
  properties = list(
    model = S7::new_property(S7::class_character, default = "Ex1"),
    nugget = S7::new_property(S7::class_logical, default = FALSE),
    unc = S7::new_property(S7::class_logical, default = TRUE),
    r_resol = S7::new_property(S7::class_numeric, default = 100),
    h_resol = S7::new_property(S7::class_numeric, default = 5),
    rs_type = S7::new_property(S7::class_character, default = "rtop"),
    hs_type = S7::new_property(S7::class_character, default = "regular"),
    cloud = S7::new_property(S7::class_logical, default = FALSE),
    amul = S7::new_property(S7::class_numeric, default = 2),
    dmul = S7::new_property(S7::class_numeric, default = 3),
    fit_method = S7::new_property(S7::class_numeric, default = 9),
    g_dist_est = S7::new_property(S7::class_logical, default = FALSE),
    g_dist_pred = S7::new_property(S7::class_logical, default = FALSE),
    var_clean = S7::new_property(S7::class_logical, default = FALSE),
    max_dist = S7::new_property(S7::class_numeric, default = Inf),
    n_max = S7::new_property(S7::class_numeric, default = 10),
    n_clus = S7::new_property(S7::class_numeric, default = 1),
    cn_areas = S7::new_property(S7::class_numeric, default = 100),
    clus_type = utop_optional(S7::class_character),
    outfile = utop_optional(S7::class_character),
    partial_overlap = S7::new_property(S7::class_logical, default = FALSE),
    wlim = S7::new_property(S7::class_numeric, default = 1.5),
    wlim_method = S7::new_property(S7::class_character, default = "all"),
    singular_solve = S7::new_property(S7::class_logical, default = FALSE),
    uk_trend_support = S7::new_property(
      S7::class_character,
      default = "centroid"
    ),
    cv = S7::new_property(S7::class_logical, default = FALSE),
    debug_level = S7::new_property(
      S7::class_numeric,
      default = if (interactive()) 1 else 0
    ),
    par_init = utop_optional()
  )
)

#' @rdname utop-s7-classes
#' @export
UtopVariogramModel <- S7::new_class(
  "UtopVariogramModel",
  package = "utop",
  properties = list(
    model = S7::new_property(S7::class_character, default = "Ex1"),
    params = S7::new_property(S7::class_numeric, default = c(1, 1, 0, 0, 1)),
    ss_err = utop_optional(S7::class_numeric),
    criterion = utop_optional()
  )
)

#' @rdname utop-s7-classes
#' @export
UtopVariogram <- S7::new_class(
  "UtopVariogram",
  package = "utop",
  properties = list(
    data = S7::new_property(S7::class_data.frame, default = data.frame())
  )
)

#' @rdname utop-s7-classes
#' @export
UtopVariogramCloud <- S7::new_class(
  "UtopVariogramCloud",
  package = "utop",
  properties = list(
    data = S7::new_property(S7::class_data.frame, default = data.frame())
  )
)

#' @rdname utop-s7-classes
#' @export
Utop <- S7::new_class(
  "Utop",
  package = "utop",
  properties = list(
    observations = utop_optional(),
    prediction_locations = utop_optional(),
    formula = utop_optional(S7::class_formula),
    params = utop_optional(UtopParams),
    variogram = utop_optional(UtopVariogram),
    variogram_cloud = utop_optional(UtopVariogramCloud),
    variogram_model = utop_optional(UtopVariogramModel),
    predictions = utop_optional(),
    simulations = utop_optional(),
    d_obs = utop_optional(),
    d_pred = utop_optional(),
    g_dist_obs = utop_optional(),
    g_dist_pred = utop_optional(),
    g_dist_pred_obs = utop_optional(),
    g_dist_bin = utop_optional(),
    d_bin = utop_optional(),
    overlap_obs = utop_optional(),
    overlap_pred_obs = utop_optional(),
    var_mat_obs = utop_optional(),
    var_mat_pred_obs = utop_optional(),
    var_mat_pred = utop_optional(),
    var_fit = utop_optional(),
    cv_info = utop_optional(),
    weight = utop_optional(),
    removed = utop_optional(),
    uk_residual = utop_optional(),
    check_vario = utop_optional(),
    method_parameters = utop_optional(S7::class_character)
  )
)

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

utop_object_names <- c(
  observations = "observations",
  prediction_locations = "predictionLocations",
  formula = "formulaString",
  predictions = "predictions",
  simulations = "simulations",
  d_obs = "dObs",
  d_pred = "dPred",
  g_dist_obs = "gDistObs",
  g_dist_pred = "gDistPred",
  g_dist_pred_obs = "gDistPredObs",
  g_dist_bin = "gDistBin",
  d_bin = "dBin",
  overlap_obs = "overlapObs",
  overlap_pred_obs = "overlapPredObs",
  var_mat_obs = "varMatObs",
  var_mat_pred_obs = "varMatPredObs",
  var_mat_pred = "varMatPred",
  var_fit = "varFit",
  cv_info = "cvInfo",
  weight = "weight",
  removed = "removed",
  uk_residual = "ukResidual",
  check_vario = "checkVario",
  method_parameters = "methodParameters"
)

utop_get <- function(x, name) {
  if (name %in% names(x)) {
    x[[name]]
  } else {
    NULL
  }
}

utop_params_from_list <- function(x) {
  if (S7::S7_inherits(x, UtopParams)) {
    return(x)
  }

  args <- list()
  for (new_name in names(utop_param_names)) {
    old_name <- unname(utop_param_names[[new_name]])
    if (old_name %in% names(x)) {
      args[[new_name]] <- x[[old_name]]
    }
  }
  do.call(UtopParams, args)
}

utop_params_to_list <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!S7::S7_inherits(x, UtopParams)) {
    return(x)
  }

  out <- list()
  for (new_name in names(utop_param_names)) {
    old_name <- unname(utop_param_names[[new_name]])
    value <- S7::prop(x, new_name)
    if (!is.null(value)) {
      out[[old_name]] <- value
    }
  }
  class(out) <- "rtopParams"
  out
}

utop_variogram_from_legacy <- function(x) {
  if (is.null(x) || S7::S7_inherits(x, UtopVariogram)) {
    return(x)
  }
  data <- x
  class(data) <- "data.frame"
  UtopVariogram(data = data)
}

utop_variogram_to_legacy <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!S7::S7_inherits(x, UtopVariogram)) {
    return(x)
  }
  data <- x@data
  class(data) <- c("rtopVariogram", "data.frame")
  data
}

utop_variogram_cloud_from_legacy <- function(x) {
  if (is.null(x) || S7::S7_inherits(x, UtopVariogramCloud)) {
    return(x)
  }
  data <- x
  class(data) <- "data.frame"
  UtopVariogramCloud(data = data)
}

utop_variogram_cloud_to_legacy <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!S7::S7_inherits(x, UtopVariogramCloud)) {
    return(x)
  }
  data <- x@data
  class(data) <- c("rtopVariogramCloud", "data.frame")
  data
}

utop_variogram_model_from_legacy <- function(x) {
  if (is.null(x) || S7::S7_inherits(x, UtopVariogramModel)) {
    return(x)
  }

  UtopVariogramModel(
    model = x$model,
    params = x$params,
    ss_err = attr(x, "SSErr", exact = TRUE),
    criterion = attr(x, "criterion", exact = TRUE)
  )
}

utop_variogram_model_to_legacy <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!S7::S7_inherits(x, UtopVariogramModel)) {
    return(x)
  }

  out <- list(model = x@model, params = x@params)
  class(out) <- "rtopVariogramModel"
  if (!is.null(x@ss_err)) {
    attr(out, "SSErr") <- x@ss_err
  }
  if (!is.null(x@criterion)) {
    attr(out, "criterion") <- x@criterion
  }
  out
}

utop_from_rtop <- function(x) {
  if (S7::S7_inherits(x, Utop)) {
    return(x)
  }

  args <- list(
    params = utop_params_from_list(utop_get(x, "params")),
    variogram = utop_variogram_from_legacy(utop_get(x, "variogram")),
    variogram_cloud = utop_variogram_cloud_from_legacy(utop_get(
      x,
      "variogramCloud"
    )),
    variogram_model = utop_variogram_model_from_legacy(utop_get(
      x,
      "variogramModel"
    ))
  )

  for (new_name in names(utop_object_names)) {
    old_name <- unname(utop_object_names[[new_name]])
    value <- utop_get(x, old_name)
    if (!is.null(value)) {
      args[[new_name]] <- value
    }
  }

  do.call(Utop, args)
}

utop_to_rtop <- function(x) {
  if (!S7::S7_inherits(x, Utop)) {
    return(x)
  }

  out <- list()
  if (!is.null(x@params)) {
    out$params <- utop_params_to_list(x@params)
  }
  if (!is.null(x@variogram)) {
    out$variogram <- utop_variogram_to_legacy(x@variogram)
  }
  if (!is.null(x@variogram_cloud)) {
    out$variogramCloud <- utop_variogram_cloud_to_legacy(x@variogram_cloud)
  }
  if (!is.null(x@variogram_model)) {
    out$variogramModel <- utop_variogram_model_to_legacy(x@variogram_model)
  }

  for (new_name in names(utop_object_names)) {
    old_name <- unname(utop_object_names[[new_name]])
    value <- S7::prop(x, new_name)
    if (!is.null(value)) {
      out[[old_name]] <- value
    }
  }

  class(out) <- "rtop"
  out
}
