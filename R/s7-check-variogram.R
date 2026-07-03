#' Create diagnostic variogram checks
#'
#' @param object A [Utop] object or [UtopVariogramModel].
#' @param ... Method-specific arguments. For [Utop] objects, supply `params`
#'   as a [UtopParams] object.
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

S7::method(utop_check_variogram, Utop) <- function(object, params = NULL, ...) {
  compute_check_variogram(object, params = params, ...)
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
