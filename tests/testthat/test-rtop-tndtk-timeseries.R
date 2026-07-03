test_that("sf variogram uses explicit support area", {
  demo <- utop_demo_data()
  observations <- demo$gauged_catchments[1:5, ]
  observations$obs <- demo$gauged_catchments$MAF[1:5]
  observations$area <- seq(10, 50, by = 10)

  vario <- utop_variogram(
    observations,
    formula = obs ~ 1,
    params = utop_params(cloud = TRUE)
  )

  expect_true(S7::S7_inherits(vario, UtopVariogramCloud))
  expect_equal(vario@data$a1, observations$area[vario@data$acl1])
  expect_equal(vario@data$a2, observations$area[vario@data$acl2])
})

test_that("stars pooled variogram averages TNDTK time steps", {
  fixture <- utop_tndtk_stars_fixture(n_obs = 6, n_pred = 2)
  observations <- fixture$observations
  obs_mat <- observations[["obs"]]

  vario <- utop_variogram(
    observations,
    formula = obs ~ 1,
    params = utop_params(cloud = TRUE)
  )
  manual_gamma <- mean((obs_mat[1, ] - obs_mat[2, ])^2 / 2)
  pair_vario <- vario@data[vario@data$acl1 == 1 & vario@data$acl2 == 2, ]
  support <- utop:::utop_stars_support(observations)

  expect_s3_class(support, "sf")
  expect_true(S7::S7_inherits(vario, UtopVariogramCloud))
  expect_equal(nrow(pair_vario), 1)
  expect_equal(pair_vario$np, length(fixture$time))
  expect_equal(pair_vario$gamma, manual_gamma, tolerance = 1e-12)
  expect_equal(pair_vario$a1, support$area[1])
  expect_equal(pair_vario$a2, support$area[2])
})

test_that("TNDTK stars time series interpolate with pooled variogram", {
  fixture <- utop_tndtk_stars_fixture(n_obs = 8, n_pred = 3)
  params <- utop_params(
    r_resol = 4,
    rs_type = "regular",
    debug_level = -1,
    nugget = FALSE,
    cloud = FALSE
  )

  rtop_obj <- utop_object(
    fixture$observations,
    fixture$prediction_locations,
    formula = obs ~ 1,
    params = params
  )
  rtop_obj <- utop_variogram(rtop_obj)

  expect_true(S7::S7_inherits(rtop_obj@variogram, UtopVariogram))
  expect_true(all(rtop_obj@variogram@data$np %% length(fixture$time) == 0))
  expect_gt(max(rtop_obj@variogram@data$np), length(fixture$time))

  rtop_obj <- utop_fit_variogram(rtop_obj, iprint = -1)
  result <- utop_krige(rtop_obj)
  predictions <- result@predictions

  expect_s3_class(predictions, "stars")
  expect_equal(unname(dim(predictions))[1], 3)
  expect_equal(unname(dim(predictions))[2], length(fixture$time))
  expect_true(all(is.finite(predictions[["var1.pred"]])))
  expect_true(all(is.finite(predictions[["var1.var"]])))
  expect_gt(stats::sd(as.vector(predictions[["var1.pred"]])), 0)
})
