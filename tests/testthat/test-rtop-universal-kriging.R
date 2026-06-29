fixtures <- utop_spatial_fixtures()

set.seed(1501)
uk_base <- createRtopObject(
  fixtures$observations,
  fixtures$prediction_locations,
  params = fixtures$params,
  formulaString = "obs ~ 1"
)
uk_base <- rtopFitVariogram(uk_base, iprint = -1)
uk_base <- rtopKrige(uk_base) # also creates varMatObs / varMatPredObs

obs_xy <- utop:::utop_centroid_coordinates(fixtures$observations)
pred_xy <- utop:::utop_centroid_coordinates(fixtures$prediction_locations)

test_that("ukTrendMatrix builds the trend basis from the RHS", {
  observations <- fixtures$observations
  observations$elev <- obs_xy[, 2] / 1e5

  f_ok <- utop:::ukTrendMatrix(obs ~ 1, observations)
  expect_equal(dim(f_ok), c(30, 1))
  expect_true(all(f_ok == 1))

  f_attr <- utop:::ukTrendMatrix(obs ~ elev, observations)
  expect_equal(dim(f_attr), c(30, 2))
  expect_equal(unname(f_attr[, 2]), observations$elev)

  f_cor <- utop:::ukTrendMatrix(obs ~ x + y, observations)
  expect_equal(unname(f_cor[, 2]), unname(obs_xy[, 1]))
  expect_equal(unname(f_cor[, 3]), unname(obs_xy[, 2]))

  expect_error(
    utop:::ukTrendMatrix(obs ~ missingVar, observations),
    "not found in data"
  )
})

test_that("block-averaged coordinate basis matches centroids for supports", {
  params_block <- getRtopParams(list(
    ukTrendSupport = "block",
    rresol = 25,
    debug.level = -1
  ))
  f_block <- utop:::ukTrendMatrix(
    obs ~ x + y,
    fixtures$observations,
    params_block
  )
  f_cent <- utop:::ukTrendMatrix(obs ~ x + y, fixtures$observations)

  # the centroid is the mean coordinate, so the block average of x and y
  # should be close to, but not identical with, the centroid evaluation
  expect_equal(unname(f_block[, 2]), unname(f_cent[, 2]), tolerance = 0.01)
  expect_equal(unname(f_block[, 3]), unname(f_cent[, 3]), tolerance = 0.01)
  expect_false(isTRUE(all.equal(unname(f_block[, 2]), unname(f_cent[, 2]))))
})

test_that("universal kriging reproduces an exact attribute trend", {
  observations <- fixtures$observations
  predictionLocations <- fixtures$prediction_locations
  observations$elev <- obs_xy[, 2] / 1e5
  predictionLocations$elev <- pred_xy[, 2] / 1e5
  observations$obs <- 2 + 3 * observations$elev

  ret <- rtopKrige(
    observations,
    predictionLocations,
    varMatObs = uk_base$varMatObs,
    varMatPredObs = uk_base$varMatPredObs,
    formulaString = obs ~ elev,
    params = list(debug.level = -1),
    wlim = Inf
  )
  expect_equal(
    ret$predictions$var1.pred,
    2 + 3 * predictionLocations$elev,
    tolerance = 1e-8
  )

  ret_cv <- rtopKrige(
    observations,
    varMatObs = uk_base$varMatObs,
    varMatPredObs = NULL,
    formulaString = obs ~ elev,
    params = list(debug.level = -1),
    cv = TRUE,
    wlim = Inf
  )
  expect_equal(
    ret_cv$predictions$var1.pred,
    2 + 3 * observations$elev,
    tolerance = 1e-8
  )
  expect_true(all(ret_cv$predictions$var1.var > 0))
})

test_that("universal kriging reproduces an exact coordinate trend", {
  observations <- fixtures$observations
  predictionLocations <- fixtures$prediction_locations
  observations$obs <- 1 + 2e-6 * obs_xy[, 1] + 3e-6 * obs_xy[, 2]

  ret <- rtopKrige(
    observations,
    predictionLocations,
    varMatObs = uk_base$varMatObs,
    varMatPredObs = uk_base$varMatPredObs,
    formulaString = obs ~ x + y,
    params = list(debug.level = -1),
    wlim = Inf
  )
  expect_equal(
    ret$predictions$var1.pred,
    unname(1 + 2e-6 * pred_xy[, 1] + 3e-6 * pred_xy[, 2]),
    tolerance = 1e-8
  )
})

test_that("block-support trend reproduces a trend in block-averaged basis", {
  observations <- fixtures$observations
  predictionLocations <- fixtures$prediction_locations
  params_block <- getRtopParams(list(
    ukTrendSupport = "block",
    rresol = 25,
    debug.level = -1
  ))
  f_obs <- utop:::ukTrendMatrix(obs ~ x, observations, params_block)
  f_pred <- utop:::ukTrendMatrix(obs ~ x, predictionLocations, params_block)
  observations$obs <- 4 + 5e-6 * f_obs[, 2]

  ret <- rtopKrige(
    observations,
    predictionLocations,
    varMatObs = uk_base$varMatObs,
    varMatPredObs = uk_base$varMatPredObs,
    formulaString = obs ~ x,
    params = list(ukTrendSupport = "block", rresol = 25, debug.level = -1),
    wlim = Inf
  )
  expect_equal(
    ret$predictions$var1.pred,
    4 + 5e-6 * unname(f_pred[, 2]),
    tolerance = 1e-8
  )
})

test_that("full rtop pipeline works with a universal kriging formula", {
  observations <- fixtures$observations
  predictionLocations <- fixtures$prediction_locations
  observations$elev <- obs_xy[, 2] / 1e5
  predictionLocations$elev <- pred_xy[, 2] / 1e5

  set.seed(1501)
  uk_obj <- createRtopObject(
    observations,
    predictionLocations,
    params = fixtures$params,
    formulaString = obs ~ elev
  )
  uk_obj <- rtopFitVariogram(uk_obj, iprint = -1)
  uk_obj <- rtopKrige(uk_obj)
  expect_true(all(is.finite(uk_obj$predictions$var1.pred)))
  expect_true(all(uk_obj$predictions$var1.var > 0))

  uk_cv <- rtopKrige(uk_obj, cv = TRUE)
  expect_true(all(is.finite(uk_cv$predictions$var1.pred)))
  expect_true(all(is.finite(uk_cv$predictions$residual)))
})

test_that("sample variogram is computed from trend residuals", {
  observations <- fixtures$observations
  observations$elev <- obs_xy[, 2] / 1e5
  trended <- observations
  trended$obs <- trended$obs + 10 * trended$elev

  vario_raw <- rtopVariogram(
    trended,
    formulaString = obs ~ 1,
    params = fixtures$params
  )
  vario_res <- rtopVariogram(
    trended,
    formulaString = obs ~ elev,
    params = fixtures$params
  )
  vario_orig <- rtopVariogram(
    observations,
    formulaString = obs ~ 1,
    params = fixtures$params
  )

  expect_lt(mean(vario_res$gamma), mean(vario_raw$gamma))
  # the residual variogram should be of the same order as the variogram of
  # the original (untrended) field
  expect_lt(mean(vario_res$gamma), 10 * mean(vario_orig$gamma))
})

test_that("spatiotemporal universal kriging reproduces an exact trend", {
  st_fixtures <- utop_stars_fixtures(n_obs = 8, n_pred = 4, n_time = 3)
  st_obs <- st_fixtures$observations
  st_pred <- st_fixtures$prediction_locations
  st_obs <- utop:::utop_stars_set_attr_matrix(
    st_obs,
    "cov",
    matrix(seq_len(8) / 2, nrow = 8, ncol = 3)
  )
  st_pred <- utop:::utop_stars_set_attr_matrix(
    st_pred,
    "cov",
    matrix(runif(4, 0, 4), nrow = 4, ncol = 3)
  )

  base_obj <- createRtopObject(
    st_obs,
    st_pred,
    formulaString = "obs ~ 1",
    params = list(
      rresol = 4,
      rstype = "regular",
      debug.level = -1,
      nugget = FALSE
    )
  )
  base_obj <- rtopFitVariogram(base_obj, iprint = -1)
  base_obj <- rtopKrige(base_obj)

  obs_cov <- utop:::utop_stars_attr_matrix(st_obs, "cov")
  pred_cov <- utop:::utop_stars_attr_matrix(st_pred, "cov")
  st_exact <- utop:::utop_stars_set_attr_matrix(st_obs, "obs", 2 + 3 * obs_cov)

  ret <- rtopKrige(
    st_exact,
    st_pred,
    varMatObs = base_obj$varMatObs,
    varMatPredObs = base_obj$varMatPredObs,
    formulaString = obs ~ cov,
    params = list(debug.level = -1),
    wlim = Inf
  )
  expect_equal(
    as.vector(ret$predictions[["var1.pred"]]),
    as.vector(2 + 3 * pred_cov),
    tolerance = 1e-8
  )

  ret_cv <- rtopKrige(
    st_exact,
    varMatObs = base_obj$varMatObs,
    formulaString = obs ~ cov,
    params = list(debug.level = -1),
    cv = TRUE,
    wlim = Inf
  )
  expect_equal(
    as.vector(ret_cv$predictions[["var1.pred"]]),
    as.vector(2 + 3 * obs_cov),
    tolerance = 1e-8
  )
})

test_that("spatiotemporal sample variogram uses trend residuals", {
  st_fixtures <- utop_stars_fixtures(n_obs = 8, n_pred = 4, n_time = 3)
  st_obs <- st_fixtures$observations
  st_obs <- utop:::utop_stars_set_attr_matrix(
    st_obs,
    "cov",
    matrix(seq_len(8) * 2, nrow = 8, ncol = 3)
  )
  cov <- utop:::utop_stars_attr_matrix(st_obs, "cov")
  obs <- utop:::utop_stars_attr_matrix(st_obs, "obs")
  st_obs <- utop:::utop_stars_set_attr_matrix(st_obs, "obs", obs + 3 * cov)

  vario_raw <- rtopVariogram(st_obs, formulaString = "obs ~ 1")
  vario_res <- rtopVariogram(st_obs, formulaString = "obs ~ cov")

  expect_s3_class(vario_res, "rtopVariogram")
  expect_lt(mean(vario_res$gamma), mean(vario_raw$gamma))
})
