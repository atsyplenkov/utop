test_that("S7 utop pipeline works for sf data", {
  rpath <- system.file("extdata", package = "utop")
  observations <- sf::st_read(rpath, "observations", quiet = TRUE)
  prediction_locations <- sf::st_read(
    rpath,
    "predictionLocations",
    quiet = TRUE
  )
  observations$obs <- observations$QSUMMER_OB / observations$AREASQKM

  obj <- utop_object(
    observations,
    prediction_locations[1:5, ],
    formula = obs ~ 1,
    params = list(g_dist_est = TRUE, g_dist_pred = TRUE, cloud = FALSE)
  )

  expect_true(S7::S7_inherits(obj, Utop))
  expect_true(S7::S7_inherits(obj@params, UtopParams))
  expect_s3_class(obj@observations, "sf")

  obj <- utop_fit_variogram(obj, iprint = -1)
  expect_true(S7::S7_inherits(obj@variogram, UtopVariogram))
  expect_true(S7::S7_inherits(obj@variogram_model, UtopVariogramModel))

  obj <- utop_krige(obj)
  expect_s3_class(obj@predictions, "sf")
  expect_equal(nrow(obj@predictions), 5)
  expect_false(anyNA(obj@predictions$var1.pred))
})

test_that("S7 variogram model constructor keeps fitting metadata", {
  model <- utop_variogram_model(nugget = 0.1)

  expect_true(S7::S7_inherits(model, UtopVariogramModel))
  expect_equal(model@model, "Ex1")
  expect_equal(model@params[3], 0.1)
})
