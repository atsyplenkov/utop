test_that("SPDF variogram uses explicit support area", {
  demo <- utop_demo_data()
  observations <- as(demo$gauged_catchments[1:5, ], "Spatial")
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

test_that("STSDF pooled variogram averages TNDTK time steps", {
  skip_if_not_installed("spacetime")

  fixture <- utop_tndtk_spacetime_fixture(n_obs = 6, n_pred = 2)
  observations <- fixture$observations
  obs <- as.data.frame(observations)
  obs$sp_index <- as.integer(as.character(obs$sp.ID))

  vario <- rtopVariogram(
    observations,
    formulaString = obs ~ 1,
    params = list(cloud = TRUE)
  )
  pair <- obs[obs$sp_index %in% c(1, 2), c("sp_index", "timeIndex", "obs")]
  pair_wide <- reshape(
    pair,
    idvar = "timeIndex",
    timevar = "sp_index",
    direction = "wide"
  )
  manual_gamma <- mean((pair_wide$obs.1 - pair_wide$obs.2)^2 / 2)
  pair_vario <- vario[vario$acl1 == 1 & vario$acl2 == 2, ]

  expect_s4_class(observations@sp, "SpatialPolygonsDataFrame")
  expect_s3_class(vario, "rtopVariogramCloud")
  expect_equal(nrow(pair_vario), 1)
  expect_equal(pair_vario$np, length(fixture$time))
  expect_equal(pair_vario$gamma, manual_gamma, tolerance = 1e-12)
  expect_equal(pair_vario$a1, observations@sp$area[1])
  expect_equal(pair_vario$a2, observations@sp$area[2])
})

test_that("TNDTK SPDF time series interpolate with pooled variogram", {
  skip_if_not_installed("spacetime")

  fixture <- utop_tndtk_spacetime_fixture(n_obs = 8, n_pred = 3)
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

  expect_s4_class(predictions, "STSDF")
  expect_equal(unname(dim(predictions))[1], 3)
  expect_equal(unname(dim(predictions))[2], length(fixture$time))
  expect_true(all(is.finite(predictions@data$var1.pred)))
  expect_true(all(is.finite(predictions@data$var1.var)))
  expect_gt(stats::sd(predictions@data$var1.pred), 0)
})
