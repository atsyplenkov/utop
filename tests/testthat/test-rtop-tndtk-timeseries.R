test_that("sf variogram uses explicit support area", {
  demo <- utop_demo_data()
  observations <- demo$gauged_catchments[1:5, ]
  observations$obs <- demo$gauged_catchments$MAF[1:5]
  observations$area <- seq(10, 50, by = 10)

  vario <- rtopVariogram(
    observations,
    formulaString = obs ~ 1,
    params = list(cloud = TRUE)
  )

  expect_s3_class(vario, "rtopVariogramCloud")
  expect_equal(vario$a1, observations$area[vario$acl1])
  expect_equal(vario$a2, observations$area[vario$acl2])
})

test_that("stars pooled variogram averages TNDTK time steps", {
  fixture <- utop_tndtk_stars_fixture(n_obs = 6, n_pred = 2)
  observations <- fixture$observations
  obs_mat <- observations[["obs"]]

  vario <- rtopVariogram(
    observations,
    formulaString = obs ~ 1,
    params = list(cloud = TRUE)
  )
  manual_gamma <- mean((obs_mat[1, ] - obs_mat[2, ])^2 / 2)
  pair_vario <- vario[vario$acl1 == 1 & vario$acl2 == 2, ]
  support <- utop:::utop_stars_support(observations)

  expect_s3_class(support, "sf")
  expect_s3_class(vario, "rtopVariogramCloud")
  expect_equal(nrow(pair_vario), 1)
  expect_equal(pair_vario$np, length(fixture$time))
  expect_equal(pair_vario$gamma, manual_gamma, tolerance = 1e-12)
  expect_equal(pair_vario$a1, support$area[1])
  expect_equal(pair_vario$a2, support$area[2])
})

test_that("TNDTK stars time series interpolate with pooled variogram", {
  fixture <- utop_tndtk_stars_fixture(n_obs = 8, n_pred = 3)
  params <- list(
    rresol = 4,
    rstype = "regular",
    debug.level = -1,
    nugget = FALSE,
    cloud = FALSE
  )

  rtop_obj <- createRtopObject(
    fixture$observations,
    fixture$prediction_locations,
    formulaString = obs ~ 1,
    params = params
  )
  rtop_obj <- rtopVariogram(rtop_obj)

  expect_s3_class(rtop_obj$variogram, "rtopVariogram")
  expect_true(all(rtop_obj$variogram$np %% length(fixture$time) == 0))
  expect_gt(max(rtop_obj$variogram$np), length(fixture$time))

  rtop_obj <- rtopFitVariogram(rtop_obj, iprint = -1)
  result <- rtopKrige(rtop_obj)
  predictions <- result$predictions

  expect_s3_class(predictions, "stars")
  expect_equal(unname(dim(predictions))[1], 3)
  expect_equal(unname(dim(predictions))[2], length(fixture$time))
  expect_true(all(is.finite(predictions[["var1.pred"]])))
  expect_true(all(is.finite(predictions[["var1.var"]])))
  expect_gt(stats::sd(as.vector(predictions[["var1.pred"]])), 0)
})
