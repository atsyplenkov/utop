fixtures <- utop_spatial_fixtures()

set.seed(1501)
rtop_fitted <- utop_object(
  fixtures$observations,
  fixtures$prediction_locations,
  params = fixtures$params,
  formula = "obs ~ 1"
)
rtop_fitted <- utop_fit_variogram(rtop_fitted, iprint = -1)

test_that("kriging returns expected prediction structures", {
  rtop_cv <- utop_krige(rtop_fitted, cv = TRUE)
  rtop_pred <- utop_krige(rtop_fitted)

  expect_true(S7::S7_inherits(rtop_fitted@variogram_model, UtopVariogramModel))

  expect_s3_class(rtop_cv@predictions, "sf")
  expect_s3_class(rtop_pred@predictions, "sf")
  expect_equal(nrow(rtop_cv@predictions), 30)
  expect_equal(nrow(rtop_pred@predictions), 2)
  expect_true(all(
    c("observed", "var1.pred", "var1.var") %in% names(rtop_cv@predictions)
  ))
  expect_true(all(c("var1.pred", "var1.var") %in% names(rtop_pred@predictions)))
})

test_that("kriging reuses the legacy semivariance path consistently", {
  rtop_cv <- utop_krige(rtop_fitted, cv = TRUE)
  rtop_pred <- utop_krige(rtop_fitted)
  varmat <- utop_var_mat(
    fixtures$observations,
    fixtures$prediction_locations,
    variogram_model = rtop_fitted@variogram_model,
    g_dist_est = TRUE,
    g_dist_pred = TRUE,
    r_resol = 25,
    h_resol = 3
  )
  rtop_reuse <- utop_krige(rtop_cv)

  expect_s3_class(rtop_reuse@predictions, "sf")
  expect_equal(nrow(rtop_reuse@predictions), 2)
  expect_true(isTRUE(all.equal(varmat@var_mat_obs, rtop_cv@var_mat_obs)))
  expect_true(isTRUE(all.equal(rtop_reuse@predictions, rtop_pred@predictions)))
})

test_that("spatial variogram updates rebuild semivariance matrices", {
  rtop_reuse <- utop_krige(utop_krige(rtop_fitted, cv = TRUE))

  rtop_updated <- utop_var_mat(rtop_reuse)
  rtop_updated <- utop_update_variogram(
    rtop_updated,
    exp = 1.5,
    action = "mult"
  )
  rtop_updated_mat <- utop_var_mat(rtop_updated)

  expect_false(isTRUE(all.equal(
    rtop_updated@var_mat_obs,
    rtop_updated_mat@var_mat_obs
  )))
  expect_true(!is.null(rtop_updated_mat@var_mat_obs))
})

test_that("spatial simulation returns expected structure", {
  set.seed(1501)
  rtop_sim <- utop_sim(
    rtop_fitted,
    nsim = 2,
    logdist = TRUE,
    debug_level = -1
  )

  expect_s3_class(rtop_sim@simulations, "data.frame")
  expect_true(all(c("sim1", "sim2") %in% names(rtop_sim@simulations)))
  expect_equal(
    nrow(rtop_sim@simulations),
    nrow(rtop_fitted@prediction_locations)
  )
})
