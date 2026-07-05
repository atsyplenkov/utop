# Regression tests for the parallel (n_clus > 1) compute paths.
#
# Parallelism in utop is a speed optimisation, not a different algorithm: the
# serial and parallel branches must produce numerically identical output. These
# tests pin that invariant by computing each result both ways and asserting
# equality. They guard against the class of bug where a parallel branch dies
# silently or produces partial/NA output.
#
# Discretisation (rs_type = "rtop") consumes from R's global RNG stream, so
# serial and parallel runs must be seeded identically before utop_var_mat to
# avoid conflating RNG drift with parallel/serial divergence.

fix <- utop_spatial_fixtures()

# A fitted variogram model shared across the var_mat tests so that parallel vs
# serial comparisons are against the same semivariance function.
fit_for_vm <- utop_fit_variogram(
  utop_object(
    fix$observations,
    fix$prediction_locations,
    params = fix$params,
    formula = "obs ~ 1"
  ),
  iprint = -1
)
vmod <- fit_for_vm@variogram_model

# Common parallel params: n_clus = 2 and a low cn_areas so the parallel gate
# (n_clus > 1 && ... > cn_areas) actually triggers. `outfile = ""` suppresses
# the per-worker "starting worker pid" console spam.
parallel_params <- function(...) {
  utop_params(n_clus = 2, cn_areas = 5, debug_level = -1, outfile = "", ...)
}

# Stop any cached cluster before and after the expression runs, so worker
# connections are closed cleanly rather than leaking unserialize warnings into
# subsequent tests. A fixed seed is set before `expr` so that RNG-dependent
# discretisation is reproducible across serial/parallel comparisons.
with_clean_cluster <- function(expr, seed = 42) {
  old <- getOption("utop.cluster")
  if (!is.null(old)) {
    parallel::stopCluster(old)
  }
  options(utop.cluster = NULL)
  set.seed(seed)
  on.exit(
    {
      cur <- getOption("utop.cluster")
      if (!is.null(cur)) {
        parallel::stopCluster(cur)
      }
      options(utop.cluster = old)
    },
    add = TRUE
  )
  force(expr)
}

# Helper: a fresh Utop with observations + prediction_locations + a fitted
# variogram model attached, sharing `vmod` so parallel/serial compare against
# the same semivariance function.
utop_for_vm <- function(params) {
  obj <- utop_object(
    fix$observations,
    fix$prediction_locations,
    params = params,
    formula = "obs ~ 1"
  )
  obj@variogram_model <- vmod
  obj
}

test_that("parallel var_mat (compute_var_mat_list path) matches serial", {
  # g_dist_* = FALSE routes through compute_var_mat_list, whose parallel branch
  # uses clusterExport + clusterApplyLB.
  serial <- with_clean_cluster(utop_var_mat(utop_for_vm(utop_params(
    g_dist_est = FALSE,
    g_dist_pred = FALSE,
    r_resol = 25,
    h_resol = 3,
    debug_level = -1
  ))))

  par <- with_clean_cluster(utop_var_mat(utop_for_vm(parallel_params(
    g_dist_est = FALSE,
    g_dist_pred = FALSE,
    r_resol = 25,
    h_resol = 3
  ))))

  expect_true(S7::S7_inherits(par, Utop))
  expect_equal(dim(par@var_mat_obs), c(30L, 30L))
  expect_true(all(is.finite(par@var_mat_obs)))
  expect_equal(par@var_mat_obs, serial@var_mat_obs)
  expect_equal(par@var_mat_pred_obs, serial@var_mat_pred_obs)
})

test_that("parallel var_mat (eval_var_mat_g_dist path) matches serial", {
  # g_dist_pred = TRUE routes through compute_var_mat_utop's g_dist branch,
  # whose parallel evaluation uses eval_var_mat_g_dist (parLapply, no
  # clusterExport). Exercises the obs-obs and obs-pred parallel eval blocks.
  serial <- with_clean_cluster(utop_var_mat(utop_for_vm(utop_params(
    g_dist_est = TRUE,
    g_dist_pred = TRUE,
    r_resol = 25,
    h_resol = 3,
    debug_level = -1
  ))))

  par <- with_clean_cluster(utop_var_mat(utop_for_vm(parallel_params(
    g_dist_est = TRUE,
    g_dist_pred = TRUE,
    r_resol = 25,
    h_resol = 3
  ))))

  expect_equal(dim(par@var_mat_obs), c(30L, 30L))
  expect_equal(dim(par@var_mat_pred_obs), c(30L, 2L))
  expect_true(all(is.finite(par@var_mat_obs)))
  expect_equal(par@var_mat_obs, serial@var_mat_obs)
  expect_equal(par@var_mat_pred_obs, serial@var_mat_pred_obs)
})

test_that("parallel discretisation (compute_disc_sf) matches serial", {
  # The parallel branch of compute_disc_sf uses clusterExport + clusterApply
  # with lfun carried as a closure. Discretisation point counts per area must
  # match exactly between serial and parallel (deterministic given same seed).
  serial <- with_clean_cluster(utop_disc(
    fix$observations,
    params = utop_params(r_resol = 100, debug_level = -1)
  ))
  par <- with_clean_cluster(utop_disc(
    fix$observations,
    params = parallel_params(r_resol = 100)
  ))

  expect_length(par, nrow(fix$observations))
  expect_length(serial, nrow(fix$observations))

  serial_lens <- vapply(
    serial,
    function(x) nrow(utop_point_coordinates(x)),
    integer(1)
  )
  par_lens <- vapply(
    par,
    function(x) nrow(utop_point_coordinates(x)),
    integer(1)
  )
  expect_identical(par_lens, serial_lens)
  expect_true(all(par_lens > 0))
})
