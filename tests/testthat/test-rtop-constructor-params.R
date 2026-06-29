test_that("constructor and parameter flow cover inference, coercion, and projection checks", {
  spatial <- utop_spatial_subset_fixtures(n_obs = 6, n_pred = 2)
  sf_fixtures <- utop_sf_subset_fixtures(n_obs = 6, n_pred = 2)

  expect_warning(
    inferred <- createRtopObject(
      spatial$observations,
      spatial$prediction_locations,
      params = spatial$params
    ),
    "formulaString missing"
  )
  expect_s3_class(inferred, "rtop")
  expect_identical(deparse(inferred$formulaString), "obs ~ 1")
  expect_true("area" %in% names(inferred$observations))
  expect_true("area" %in% names(inferred$predictionLocations))

  coerced <- createRtopObject(
    spatial$observations,
    spatial$prediction_locations,
    formulaString = "obs ~ 1",
    params = spatial$params
  )
  expect_s3_class(coerced$formulaString, "formula")
  expect_identical(deparse(coerced$formulaString), "obs ~ 1")

  updated <- createRtopObject(
    coerced,
    params = list(gDist = FALSE, model = "Ex1")
  )
  expect_s3_class(updated, "rtop")
  expect_false(isTRUE(updated$params$gDistEst))
  expect_false(isTRUE(updated$params$gDistPred))

  gdist_true <- createRtopObject(
    spatial$observations,
    spatial$prediction_locations,
    formulaString = "obs ~ 1",
    params = list(gDist = TRUE)
  )
  expect_true(gdist_true$params$gDistEst)
  expect_true(gdist_true$params$gDistPred)

  gdist_false <- createRtopObject(
    spatial$observations,
    spatial$prediction_locations,
    formulaString = "obs ~ 1",
    params = list(gDist = FALSE)
  )
  expect_false(isTRUE(gdist_false$params$gDistEst))
  expect_false(isTRUE(gdist_false$params$gDistPred))

  expect_error(utop:::getRtopParams(list(), geoDist = TRUE))

  supported_model <- createRtopObject(
    spatial$observations,
    spatial$prediction_locations,
    formulaString = "obs ~ 1",
    params = list(model = "Ex1")
  )
  expect_identical(supported_model$params$model, "Ex1")

  expect_error(
    createRtopObject(
      spatial$observations,
      spatial$prediction_locations,
      formulaString = "obs ~ 1",
      params = list(model = "Bogus")
    ),
    "not implemented"
  )

  sf_obj <- createRtopObject(
    sf_fixtures$observations,
    sf_fixtures$prediction_locations,
    formulaString = "obs ~ 1"
  )
  expect_true("area" %in% names(sf_obj$observations))
  expect_true("area" %in% names(sf_obj$predictionLocations))

  sf_obs <- sf_fixtures$observations
  sf_pred <- sf_fixtures$prediction_locations

  obs_crs <- sf::st_crs(sf_obs)
  mismatch_crs <- if (!is.na(obs_crs$epsg) && obs_crs$epsg == 4326) {
    3857
  } else {
    4326
  }

  expect_error(
    createRtopObject(
      sf_obs,
      sf::st_transform(sf_pred, mismatch_crs),
      formulaString = "obs ~ 1"
    ),
    "different projections"
  )

  expect_error(
    createRtopObject(
      sf_obs,
      sf::st_set_crs(sf_pred, NA),
      formulaString = "obs ~ 1"
    ),
    "only one of observations and predictionLocations have projection"
  )
})
