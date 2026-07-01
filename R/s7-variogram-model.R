#' @noRd
build_variogram_model <- function(
  model = "Ex1",
  sill = NULL,
  range = NULL,
  exp = NULL,
  nugget = NULL,
  exp0 = NULL,
  observations = NULL,
  formula = obs ~ 1
) {
  if (tolower(model) != "ex1") {
    stop(paste("model", model, "not implemented"), call. = FALSE)
  }
  model <- "Ex1"
  if (!is.null(observations)) {
    par_init <- find_par_init(formula, observations, model)$par0
    if (is.null(sill)) {
      sill <- par_init[1]
    }
    if (is.null(range)) {
      range <- par_init[2]
    }
    if (is.null(nugget)) {
      nugget <- par_init[3]
    }
    if (is.null(exp)) {
      exp <- par_init[4]
    }
    if (is.null(exp0)) {
      exp0 <- par_init[5]
    }
  } else {
    if (is.null(sill)) {
      sill <- 1
    }
    if (is.null(range)) {
      range <- 1
    }
    if (is.null(exp)) {
      exp <- 0
    }
    if (is.null(nugget)) {
      nugget <- 0
    }
    if (is.null(exp0)) {
      exp0 <- 1
    }
  }

  UtopVariogramModel(
    model = model,
    params = c(sill, range, nugget, exp, exp0)
  )
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
  build_variogram_model(
    model = model,
    sill = sill,
    range = range,
    exp = exp,
    nugget = nugget,
    exp0 = exp0,
    observations = observations,
    formula = formula
  )
}

#' @noRd
update_variogram_model <- function(
  model,
  action = "mult",
  ...,
  check_vario = FALSE,
  sample_variogram = NULL,
  observations = NULL
) {
  dots <- list(...)
  if (model@model != "Ex1") {
    return(model)
  }

  if ("sill" %in% names(dots)) {
    value <- dots$sill
    if (action %in% c("mult", "add")) {
      model@params[1] <- model@params[1] * value
    } else if (action == "replace") {
      model@params[1] <- value
    }
  }
  if ("range" %in% names(dots)) {
    value <- dots$range
    if (action %in% c("mult", "add")) {
      model@params[2] <- model@params[2] * value
    } else if (action == "replace") {
      model@params[2] <- value
    }
  }
  if ("nugget" %in% names(dots)) {
    value <- dots$nugget
    if (action %in% c("mult", "add")) {
      model@params[3] <- model@params[3] * value
    } else if (action == "replace") {
      model@params[3] <- value
    }
  }
  if ("exp" %in% names(dots)) {
    value <- dots$exp
    if (action %in% c("mult", "add")) {
      model@params[4] <- model@params[4] * value
    } else if (action == "replace") {
      model@params[4] <- value
    }
  }
  if ("exp0" %in% names(dots)) {
    value <- dots$exp0
    if (action %in% c("mult", "add")) {
      model@params[5] <- model@params[5] * value
    } else if (action == "replace") {
      model@params[5] <- value
    }
  }

  if (isTRUE(check_vario)) {
    utop_check_variogram(
      model,
      sample_variogram = sample_variogram,
      observations = observations
    )
  }

  model
}

#' Update a variogram model
#'
#' @param object A [Utop] or [UtopVariogramModel] object.
#' @param ... Variogram parameter updates.
#'
#' @return An updated object.
#' @export
utop_update_variogram <- S7::new_generic(
  name = "utop_update_variogram",
  dispatch_args = "object",
  fun = function(object, ...) {
    S7::S7_dispatch()
  }
)

S7::method(utop_update_variogram, Utop) <- function(
  object,
  action = "mult",
  ...,
  check_vario = FALSE
) {
  sample_variogram <- object@variogram
  if (is.null(sample_variogram)) {
    sample_variogram <- object@variogram_cloud
  }
  object@variogram_model <- update_variogram_model(
    object@variogram_model,
    action = action,
    ...,
    check_vario = check_vario,
    sample_variogram = sample_variogram,
    observations = object@observations
  )
  object
}

S7::method(utop_update_variogram, UtopVariogramModel) <- function(
  object,
  action = "mult",
  ...,
  check_vario = FALSE,
  sample_variogram = NULL,
  observations = NULL
) {
  update_variogram_model(
    object,
    action = action,
    ...,
    check_vario = check_vario,
    sample_variogram = sample_variogram,
    observations = observations
  )
}
