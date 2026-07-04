# Variogram model registry.
#
# Single source of truth for the variogram models that utop supports. The
# Fortran semivariance routines in src/vred.f identify models by integer code
# (see imodel() below); the per-model parameter count, optimiser starting-value
# rows, and any implicit constraint on those parameters all live here so that
# adding or retiring a model is a one-line change.

# Default lower/upper/initial rows for the five-parameter Ex1 family.
# Shorter models (Exp, Sph, Gau, Sp1, Fra) take a leading slice of these rows,
# with Fra overriding the second (range) row. par0 (geometric mean) is filled
# in lazily by find_par_init_default() because it depends only on parl/paru.
utop_par_init_template <- function() {
  data.frame(
    parl = c(1e-06, 1e-02, 1.0e-01, 1e-5, 1e-01),
    paru = c(5.0e+02, 1.0e7, 1.0e+07, 1.5, 1.7)
  )
}

# Implicit-constraint predicate used by the SCE-UA optimiser for models whose
# parameters are not independent. Ex1 requires 2*exp + exp0 <= 1; other models
# have no such constraint.
utop_ex1_implicit <- function(pars) (2 * pars[4] + pars[5]) > 1

# User-facing models. `imodel_code` matches the Fortran case numbers in
# src/vred.f. Note that Gho (code 5) is intentionally NOT a user-facing model:
# it is an internal sentinel used by the geostatistical-distance path
# (compute_g_dist_list) and must stay out of this registry.
#
# n_pars        - number of optimised parameters (leading slice of par_init)
# par_init_rows - optional override rows (named list, e.g. list(`2` = c(...)))
# implicit      - optional constraint predicate consumed by the optimiser
utop_models <- function() {
  list(
    Exp = list(imodel_code = 1L, n_pars = 3L),
    Ex1 = list(imodel_code = 2L, n_pars = 5L, implicit = utop_ex1_implicit),
    Gau = list(imodel_code = 3L, n_pars = 3L),
    Sph = list(imodel_code = 6L, n_pars = 3L),
    Sp1 = list(imodel_code = 7L, n_pars = 4L),
    Fra = list(
      imodel_code = 8L,
      n_pars = 3L,
      par_init_rows = list(`2` = c(1e-6, 2, 0.01))
    )
  )
}

#' @noRd
utop_model <- function(name) {
  models <- utop_models()
  out <- models[[name]]
  if (is.null(out)) {
    stop(paste("model", name, "not implemented"), call. = FALSE)
  }
  out
}

#' @noRd
utop_user_models <- function() {
  names(utop_models())
}

# Fortran integer code for a model name. Includes Gho as a non-user sentinel
# (see utop_models() for the rationale); the numbers must match the
# Fortran routines in src/vred.f.
#' @noRd
imodel <- function(model) {
  if (model == "Gho") {
    return(5L)
  }
  utop_model(model)$imodel_code
}

# Implicit-constraint predicate for a model, or NULL if it has none. Used by
# the SCE-UA optimiser in fit_variogram_binned() / fit_variogram_cloud().
#' @noRd
utop_model_implicit <- function(name) {
  implicit <- utop_model(name)$implicit
  if (is.null(implicit)) NULL else implicit
}

# Default par_init (lower/upper/initial) for a model, without any
# data-driven tuning. Used when no observations are available. Per-row
# overrides (e.g. Fra's range row) are applied AFTER par0 is computed so that
# an override can pin par0 explicitly, matching the legacy semantics.
#' @noRd
find_par_init_default <- function(model) {
  spec <- utop_model(model)
  par_init <- utop_par_init_template()
  par_init$par0 <- 10**(0.5 * (log10(par_init$paru) + log10(par_init$parl)))
  par_init <- utop_apply_par_init_overrides(par_init, spec)
  par_init[seq_len(spec$n_pars), ]
}

# Apply per-row overrides (e.g. Fra's range row) to a par_init template.
# `spec` is one entry from utop_models().
#' @noRd
utop_apply_par_init_overrides <- function(par_init, spec) {
  overrides <- spec$par_init_rows
  if (is.null(overrides) || length(overrides) == 0L) {
    return(par_init)
  }
  for (row in names(overrides)) {
    par_init[as.integer(row), ] <- overrides[[row]]
  }
  par_init
}
