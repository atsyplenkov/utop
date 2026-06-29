#' Create a utop object
#'
#' `utop_object()` creates the formal S7 workflow object used by utop's
#' snake_case API.
#'
#' @param observations `sf` polygons or a vector `stars` data cube with
#'   observations.
#' @param prediction_locations `sf` polygons or a vector `stars` data cube with
#'   prediction locations.
#' @param formula Formula describing the dependent variable and trend.
#' @param params A list or [UtopParams] object with model parameters.
#' @param overlap_obs,overlap_pred_obs Optional overlap matrices.
#' @param ... Additional arguments passed to the internal constructor.
#'
#' @return A [Utop] object.
#' @export
utop_object <- function(
  observations,
  prediction_locations,
  formula,
  params = list(),
  overlap_obs = NULL,
  overlap_pred_obs = NULL,
  ...
) {
  args <- list(
    observations = observations,
    params = utop_params_to_list(params),
    ...
  )
  if (!missing(prediction_locations)) {
    args$predictionLocations <- prediction_locations
  }
  if (!missing(formula)) {
    args$formulaString <- formula
  }
  if (!is.null(overlap_obs)) {
    args$overlapObs <- overlap_obs
  }
  if (!is.null(overlap_pred_obs)) {
    args$overlapPredObs <- overlap_pred_obs
  }

  utop_from_rtop(do.call(createRtopObject, args))
}

#' Create utop parameters
#'
#' @param params Existing parameters as a list or [UtopParams] object.
#' @param new_params Parameter values to update.
#' @param observations Optional observations used to infer starting values.
#' @param formula Optional model formula.
#' @param ... Individual parameter values.
#'
#' @return A [UtopParams] object.
#' @export
utop_params <- function(
  params = list(),
  new_params = list(),
  observations = NULL,
  formula = NULL,
  ...
) {
  args <- list(
    params = utop_params_to_list(params),
    newPar = utop_params_to_list(new_params),
    ...
  )
  if (!is.null(observations)) {
    args$observations <- observations
  }
  if (!is.null(formula)) {
    args$formulaString <- formula
  }

  utop_params_from_list(do.call(getRtopParams, args))
}

#' Create or update a utop variogram model
#'
#' @param model Variogram model name.
#' @param sill,range,exp,nugget,exp0 Variogram model parameters.
#' @param observations Optional observations used to infer missing parameters.
#' @param formula Model formula.
#'
#' @return A [UtopVariogramModel] object.
#' @export
utop_variogram_model <- function(
  model = "Ex1",
  sill = NULL,
  range = NULL,
  exp = NULL,
  nugget = NULL,
  exp0 = NULL,
  observations = NULL,
  formula = obs ~ 1
) {
  utop_variogram_model_from_legacy(rtopVariogramModel(
    model = model,
    sill = sill,
    range = range,
    exp = exp,
    nugget = nugget,
    exp0 = exp0,
    observations = observations,
    formulaString = formula
  ))
}

#' Create a sample variogram
#'
#' @param object A [Utop] object.
#' @param ... Arguments passed to the variogram implementation.
#'
#' @return A [Utop] object with a variogram attached.
#' @export
utop_variogram <- S7::new_generic("utop_variogram", "object")

S7::method(utop_variogram, Utop) <- function(object, ...) {
  utop_from_rtop(rtopVariogram(utop_to_rtop(object), ...))
}

#' Fit a variogram model
#'
#' @param object A [Utop] object.
#' @param ... Arguments passed to the variogram fitting implementation.
#'
#' @return A [Utop] object with a fitted variogram model attached.
#' @export
utop_fit_variogram <- S7::new_generic("utop_fit_variogram", "object")

S7::method(utop_fit_variogram, Utop) <- function(object, ...) {
  utop_from_rtop(rtopFitVariogram(utop_to_rtop(object), ...))
}

#' Krige with a utop object
#'
#' @param object A fitted [Utop] object.
#' @param ... Arguments passed to the kriging implementation.
#'
#' @return A [Utop] object with predictions attached.
#' @export
utop_krige <- S7::new_generic("utop_krige", "object")

S7::method(utop_krige, Utop) <- function(object, ...) {
  utop_from_rtop(rtopKrige(utop_to_rtop(object), ...))
}
