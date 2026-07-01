#' @noRd
utop_check_variogram_state <- function(object) {
  state <- list(
    observations = object@observations,
    formulaString = object@formula,
    params = object@params
  )
  if (!is.null(object@prediction_locations)) {
    state$predictionLocations <- object@prediction_locations
  }
  if (!is.null(object@variogram_model)) {
    state$variogramModel <- object@variogram_model
  }
  if (!is.null(object@variogram)) {
    state$variogram <- utop_tag_variogram_class(
      utop_variogram_data(object@variogram),
      cloud = FALSE
    )
  }
  if (!is.null(object@variogram_cloud)) {
    state$variogramCloud <- utop_tag_variogram_class(
      utop_variogram_data(object@variogram_cloud),
      cloud = TRUE
    )
  }
  if (!is.null(object@var_fit)) {
    state$varFit <- object@var_fit
  }
  if (!is.null(object@d_obs)) {
    state$dObs <- object@d_obs
  }
  if (!is.null(object@g_dist_obs)) {
    state$gDistObs <- object@g_dist_obs
  }
  state
}

#' @noRd
compute_check_variogram_utop <- function(object, params = list(), ...) {
  object@params <- utop_params(object@params, new_params = params, ...)
  state <- utop_check_variogram_state(object)
  state <- compute_check_variogram(state, params = list(), ...)
  if ("checkVario" %in% names(state)) {
    object@check_vario <- state$checkVario
  }
  object
}

#' Create diagnostic variogram checks
#'
#' @param object A [Utop] object or [UtopVariogramModel].
#' @param ... Method-specific arguments. For [Utop] objects, `params` updates
#'   parameters.
#'
#' @return An updated diagnostic object.
#' @export
utop_check_variogram <- S7::new_generic(
  name = "utop_check_variogram",
  dispatch_args = "object",
  fun = function(object, ...) {
    S7::S7_dispatch()
  }
)

S7::method(utop_check_variogram, Utop) <- function(
  object,
  params = list(),
  ...
) {
  compute_check_variogram_utop(object, params = params, ...)
}

S7::method(utop_check_variogram, UtopVariogramModel) <- function(
  object,
  sample_variogram = NULL,
  observations = NULL,
  ...
) {
  sample_variogram_data <- NULL
  cloud <- FALSE
  if (S7::S7_inherits(sample_variogram, UtopVariogram)) {
    sample_variogram_data <- utop_variogram_data(sample_variogram)
  } else if (S7::S7_inherits(sample_variogram, UtopVariogramCloud)) {
    sample_variogram_data <- utop_variogram_data(sample_variogram)
    cloud <- TRUE
  } else if (!is.null(sample_variogram)) {
    sample_variogram_data <- sample_variogram
    cloud <- variogram_is_cloud(sample_variogram)
  }
  if (!is.null(sample_variogram_data)) {
    sample_variogram_data <- utop_tag_variogram_class(
      sample_variogram_data,
      cloud = cloud
    )
  }
  compute_check_variogram_model(
    object,
    sampleVariogram = sample_variogram_data,
    observations = observations,
    ...
  )
}
