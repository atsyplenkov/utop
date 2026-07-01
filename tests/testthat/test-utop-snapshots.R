test_that("spatial workflow numerical outputs are stable", {
  fixtures <- utop_spatial_fixtures()
  set.seed(1501)
  obj <- utop_object(
    fixtures$observations,
    fixtures$prediction_locations,
    params = fixtures$params,
    formula = "obs ~ 1"
  )
  obj <- utop_fit_variogram(obj, iprint = -1)

  cv <- utop_krige(obj, cv = TRUE)
  pred <- utop_krige(obj)

  expect_snapshot_numeric(
    cor(cv@predictions$observed, cv@predictions$var1.pred),
    variant = "spatial-cv-cor"
  )
  expect_snapshot_numeric(
    pred@predictions$var1.pred,
    variant = "spatial-pred"
  )
  expect_snapshot_numeric(
    obj@variogram@data$gamma,
    variant = "spatial-vario-gamma"
  )
  expect_snapshot_numeric(
    obj@variogram_model@params,
    variant = "spatial-vario-model-params"
  )
})

test_that("sf subset workflow numerical outputs are stable", {
  fixtures <- utop_sf_subset_fixtures(n_obs = 12, n_pred = 5)
  set.seed(1)
  obj <- utop_object(
    fixtures$observations,
    fixtures$prediction_locations,
    params = fixtures$params,
    formula = "obs ~1"
  )
  obj <- utop_fit_variogram(obj, iprint = -1)
  pred <- utop_krige(obj)
  cv <- utop_krige(obj, cv = TRUE)

  expect_snapshot_numeric(
    pred@predictions$var1.pred,
    variant = "sf-pred"
  )
  expect_snapshot_numeric(
    cor(cv@predictions$observed, cv@predictions$var1.pred),
    variant = "sf-cv-cor"
  )
})

test_that("spatial simulation outputs are stable with fixed seed", {
  fixtures <- utop_spatial_fixtures()
  set.seed(1501)
  obj <- utop_object(
    fixtures$observations,
    fixtures$prediction_locations,
    params = fixtures$params,
    formula = "obs ~ 1"
  )
  obj <- utop_fit_variogram(obj, iprint = -1)

  set.seed(1501)
  sim_2 <- utop_sim(obj, nsim = 2, logdist = TRUE, debug_level = -1)
  sim_input <- obj
  sim_input@prediction_locations <- sim_input@observations
  sim_input@observations$unc <- var(sim_2@observations$obs) *
    min(sim_2@observations$area) /
    sim_2@observations$area
  sim_input@prediction_locations$replaceNumber <- seq_len(nrow(
    sim_input@prediction_locations
  ))
  sim_3 <- utop_sim(
    sim_input,
    nsim = 3,
    replace = TRUE,
    debug_level = -1
  )

  expect_snapshot_numeric(
    c(
      sim_2@simulations$sim1[1],
      sim_2@simulations$sim2[1],
      sim_2@simulations$sim2[2]
    ),
    variant = "spatial-sim-2"
  )
  expect_snapshot_numeric(
    c(
      sim_3@simulations$sim1[1],
      sim_3@simulations$sim2[1],
      sim_3@simulations$sim3[14]
    ),
    variant = "spatial-sim-3"
  )
})

test_that("covariance matrix blocks are stable", {
  fixtures <- utop_spatial_subset_fixtures(n_obs = 8, n_pred = 2)
  params <- modifyList(
    fixtures$params,
    list(nugget = FALSE, model = "Ex1", g_dist_pred = FALSE)
  )
  obj <- utop_object(
    fixtures$observations,
    fixtures$prediction_locations,
    params = params,
    formula = "obs ~ 1"
  )
  obj <- utop_fit_variogram(obj, iprint = -1)
  obj <- utop_var_mat(obj)

  expect_snapshot_numeric(
    diag(obj@var_mat_obs),
    variant = "var-mat-obs-diag"
  )
  expect_snapshot_numeric(
    obj@var_mat_pred_obs[, 1],
    variant = "var-mat-pred-obs-col1"
  )
})

test_that("stars workflow numerical outputs are stable", {
  fixtures <- utop_stars_fixtures(n_obs = 8, n_pred = 4, n_time = 3)
  obj <- utop_object(
    fixtures$observations,
    fixtures$prediction_locations,
    formula = "obs ~ 1",
    params = list(
      r_resol = 4,
      rs_type = "regular",
      debug_level = -1,
      nugget = FALSE
    )
  )
  obj <- utop_fit_variogram(obj, iprint = -1)
  pred <- utop_krige(obj)

  expect_snapshot_numeric(
    as.vector(pred@predictions[["var1.pred"]])[seq_len(4)],
    variant = "stars-pred-head"
  )
})

test_that("universal kriging workflow numerical outputs are stable", {
  fixtures <- utop_spatial_fixtures()
  obs_xy <- utop:::utop_centroid_coordinates(fixtures$observations)
  pred_xy <- utop:::utop_centroid_coordinates(fixtures$prediction_locations)
  observations <- fixtures$observations
  prediction_locations <- fixtures$prediction_locations
  observations$elev <- obs_xy[, 2] / 1e5
  prediction_locations$elev <- pred_xy[, 2] / 1e5

  set.seed(1501)
  obj <- utop_object(
    observations,
    prediction_locations,
    params = fixtures$params,
    formula = "obs ~ elev"
  )
  obj <- utop_fit_variogram(obj, iprint = -1)
  pred <- utop_krige(obj)
  cv <- utop_krige(obj, cv = TRUE)

  expect_snapshot_numeric(
    pred@predictions$var1.pred,
    variant = "uk-pred"
  )
  expect_snapshot_numeric(
    cor(cv@predictions$observed, cv@predictions$var1.pred),
    variant = "uk-cv-cor"
  )
  expect_snapshot_numeric(
    obj@variogram@data$gamma,
    variant = "uk-vario-gamma"
  )
})
