#' @noRd
find_par_init_default <- function(model) {
  par_init <- data.frame(
    parl = c(1e-06, 1e-02, 1.0e-01, 1e-5, 1e-01),
    paru = c(5.0e+02, 1.0e7, 1.0e+07, 1.5, 1.7)
  )
  par_init$par0 <- 10**(0.5 * (log10(par_init$paru) + log10(par_init$parl)))

  if (model %in% c("Exp", "Sph", "Gau")) {
    par_init <- par_init[1:3, ]
  } else if (model == "Sp1") {
    par_init <- par_init[1:4, ]
  } else if (model == "Ex1") {
    par_init <- par_init
  } else if (model == "Fra") {
    par_init[2, ] <- c(1e-6, 2, 0.01)
    par_init <- par_init[1:3, ]
  } else {
    stop(paste("model", model, "not implemented"), call. = FALSE)
  }
  par_init
}

#' @noRd
find_par_init <- function(formula, observations, model) {
  observations <- utop_add_area(observations)
  if (inherits(observations, "stars")) {
    ntime <- utop_stars_ntime(observations)
    if (ntime > 20) {
      observations <- utop_stars_slice_time(
        observations,
        sample(seq_len(ntime), 20)
      )
    }
    vario <- utop_variogram_data(utop_variogram(
      observations,
      formula = formula
    ))
    a_obs <- utop_area(observations)
  } else {
    formula_use <- formula
    if (hasUkTrend(formula)) {
      observations$ukResidual <- ukResiduals(formula, observations)
      formula_use <- ukResidual ~ 1
    }
    vario <- gstat::variogram(formula_use, observations)
    a_obs <- observations$area
  }
  par_init <- data.frame(parl = c(1:5), paru = 1, par0 = 1)
  par_init[1, 1] <- min(vario$gamma) / 10
  par_init[1, 2] <- max(vario$gamma) * 500
  par_init[2, 1] <- sqrt(min(a_obs)) / 4
  par_init[2, 2] <- max(vario$dist) * 10
  minla <- min(a_obs)
  maxla <- (max(a_obs)^1.5) * max(vario$gamma)
  par_init[3, 1] <- min(vario$gamma) * minla / 100
  par_init[3, 2] <- max(vario$gamma) * maxla
  par_init[4, 1] <- 1e-5
  par_init[4, 2] <- 1.5
  par_init[5, 1] <- 0.1
  par_init[5, 2] <- 1.7
  if (model == "Ex1") {
    par_init[4, 2] <- 1
    par_init[5, 2] <- 1
  }

  par_init[, 3] <- sqrt(par_init[, 1] * par_init[, 2])
  if (model %in% c("Exp", "Sph", "Gau")) {
    par_init <- par_init[1:3, ]
  } else if (model == "Sp1") {
    par_init <- par_init[1:4, ]
  } else if (model == "Ex1") {
    par_init <- par_init
  } else if (model == "Fra") {
    par_init[2, ] <- c(1e-6, 2, 0.01)
    par_init <- par_init[1:3, ]
  } else {
    stop(paste("model", model, "not implemented"), call. = FALSE)
  }
  par_init
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
  params = NULL,
  new_params = list(),
  observations = NULL,
  formula = NULL,
  ...
) {
  values <- utop_params_apply_aliases(list(...))

  if (S7::S7_inherits(params, UtopParams)) {
    out <- params
  } else {
    out <- UtopParams()
    if (is.list(params) && length(params) > 0L) {
      out <- update_utop_params(out, utop_params_apply_aliases(params))
    }
  }

  out <- update_utop_params(out, utop_params_apply_aliases(new_params))
  out <- update_utop_params(out, values)

  if (is.null(out@par_init)) {
    if (!is.null(observations)) {
      if (is.null(formula)) {
        formula <- as.formula(utop_default_formula(observations))
      }
      out@par_init <- find_par_init(
        formula = formula,
        observations = observations,
        model = out@model
      )
    } else {
      out@par_init <- find_par_init_default(out@model)
    }
  }

  out
}
