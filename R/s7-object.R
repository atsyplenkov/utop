#' Create a utop object
#'
#' `utop_object()` creates the formal S7 workflow object used by utop's
#' snake_case API.
#'
#' @param observations `sf` polygons, a vector `stars` data cube, or an
#'   existing [Utop] object to update.
#' @param prediction_locations `sf` polygons or a vector `stars` data cube with
#'   prediction locations.
#' @param formula Formula describing the dependent variable and trend.
#' @param params A list or [UtopParams] object with model parameters.
#' @param overlap_obs,overlap_pred_obs Optional overlap matrices.
#' @param ... Additional parameter values passed to [utop_params()].
#'
#' @return A [Utop] object.
#' @export
utop_object <- function(
  observations,
  prediction_locations = NULL,
  formula = NULL,
  params = list(),
  overlap_obs = NULL,
  overlap_pred_obs = NULL,
  ...
) {
  if (S7::S7_inherits(observations, Utop)) {
    object <- observations
    object@params <- utop_params(
      params = object@params,
      new_params = c(params, list(...)),
      observations = object@observations,
      formula = if (is.null(formula)) object@formula else formula
    )
    if (!is.null(formula)) {
      object@formula <- as.formula(formula)
    }
    return(object)
  }

  if (missing(observations)) {
    stop("observations are missing", call. = FALSE)
  }

  observations <- utop_add_area(observations)

  if (!is.null(prediction_locations)) {
    prediction_locations <- utop_add_area(prediction_locations)
    check_matching_crs(observations, prediction_locations)
  }

  missing_formula <- is.null(formula)
  if (missing_formula) {
    formula <- as.formula(utop_default_formula(observations))
    warning("formula missing, using ", deparse(formula), call. = FALSE)
  }
  if (!inherits(formula, "formula")) {
    formula <- as.formula(formula)
  }

  params <- utop_params(
    params = params,
    observations = observations,
    formula = formula,
    ...
  )

  if (params@nugget && is.null(overlap_obs)) {
    overlap_obs <- find_overlap(
      observations,
      observations,
      partial_overlap = params@partial_overlap
    )
    if (!is.null(prediction_locations) && is.null(overlap_pred_obs)) {
      overlap_pred_obs <- find_overlap(
        observations,
        prediction_locations,
        partial_overlap = params@partial_overlap
      )
    }
  }

  Utop(
    observations = observations,
    prediction_locations = prediction_locations,
    formula = formula,
    params = params,
    overlap_obs = overlap_obs,
    overlap_pred_obs = overlap_pred_obs
  )
}
