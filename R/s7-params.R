# Note: find_par_init_default() is defined in R/models.R, where it shares the
# per-model par_init template and override logic with the data-driven variant
# below.

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

  # Data-driven lower/upper bounds for the five Ex1-family parameters.
  par_init <- data.frame(parl = c(1:5), paru = 1, par0 = 1)
  par_init[1, 1] <- min(vario$gamma) / 10
  par_init[1, 2] <- max(vario$gamma) * 500
  par_init[2, 1] <- sqrt(min(a_obs)) / 4
  par_init[2, 2] <- max(vario$dist) * 10
  minla <- min(a_obs)
  maxla <- (max(a_obs)^1.5) * max(vario$gamma)
  par_init[3, 1] <- min(vario$gamma) * minla / 100
  par_init[3, 2] <- max(vario$gamma) * maxla
  par_init[4, ] <- c(1e-5, 1.5, NA)
  par_init[5, ] <- c(0.1, 1.7, NA)
  if (model == "Ex1") {
    par_init[4, 2] <- 1
    par_init[5, 2] <- 1
  }
  par_init$par0 <- sqrt(par_init$parl * par_init$paru)

  spec <- utop_model(model)
  par_init <- utop_apply_par_init_overrides(par_init, spec)
  par_init[seq_len(spec$n_pars), ]
}

#' Create utop parameters
#'
#' @param ... Canonical utop parameter values supplied as named arguments.
#' @param observations Optional observations used to infer starting values.
#' @param formula Optional model formula.
#'
#' @return A [UtopParams] object.
#' @export
utop_params <- function(..., observations = NULL, formula = NULL) {
  values <- utop_require_named_values(list(...))
  out <- UtopParams()
  out <- update_utop_params(out, values)
  utop_complete_params(out, observations = observations, formula = formula)
}
