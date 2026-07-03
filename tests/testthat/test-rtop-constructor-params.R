test_that("constructor and strict params flow cover inference and CRS checks", {
  spatial <- utop_spatial_subset_fixtures(n_obs = 6, n_pred = 2)
  sf_fixtures <- utop_sf_subset_fixtures(n_obs = 6, n_pred = 2)

  expect_warning(
    inferred <- utop_object(
      spatial$observations,
      spatial$prediction_locations,
      params = spatial$params
    ),
    "formula missing"
  )
  expect_true(S7::S7_inherits(inferred, Utop))
  expect_identical(deparse(inferred@formula), "obs ~ 1")
  expect_true("area" %in% names(inferred@observations))
  expect_true("area" %in% names(inferred@prediction_locations))

  coerced <- utop_object(
    spatial$observations,
    spatial$prediction_locations,
    formula = "obs ~ 1",
    params = spatial$params
  )
  expect_s3_class(coerced@formula, "formula")
  expect_identical(deparse(coerced@formula), "obs ~ 1")

  updated <- utop_object(
    coerced,
    params = utop_test_params(
      g_dist_est = FALSE,
      g_dist_pred = FALSE,
      model = "Ex1"
    )
  )
  expect_true(S7::S7_inherits(updated, Utop))
  expect_false(isTRUE(updated@params@g_dist_est))
  expect_false(isTRUE(updated@params@g_dist_pred))

  expect_error(
    utop_object(
      spatial$observations,
      spatial$prediction_locations,
      formula = "obs ~ 1",
      params = list(g_dist_est = TRUE)
    ),
    "UtopParams"
  )
  expect_error(utop_params(gDist = TRUE), "unknown utop parameter")
  expect_error(utop_params(list()), "only accepts named arguments")

  supported_model <- utop_object(
    spatial$observations,
    spatial$prediction_locations,
    formula = "obs ~ 1",
    params = utop_params(model = "Ex1")
  )
  expect_identical(supported_model@params@model, "Ex1")
  expect_error(utop_params(model = "Bogus"), "not implemented")

  params_no_obs <- utop_test_params()
  expect_true(utop_is_default_par_init(params_no_obs@par_init, "Ex1"))
  obj <- utop_object(
    spatial$observations,
    spatial$prediction_locations,
    formula = "obs ~ 1",
    params = params_no_obs
  )
  expect_false(utop_is_default_par_init(obj@params@par_init, "Ex1"))
  expect_false(identical(obj@params@par_init, params_no_obs@par_init))

  sf_obj <- utop_object(
    sf_fixtures$observations,
    sf_fixtures$prediction_locations,
    formula = "obs ~ 1"
  )
  expect_true("area" %in% names(sf_obj@observations))
  expect_true("area" %in% names(sf_obj@prediction_locations))

  sf_obs <- sf_fixtures$observations
  sf_pred <- sf_fixtures$prediction_locations

  obs_crs <- sf::st_crs(sf_obs)
  mismatch_crs <- if (!is.na(obs_crs$epsg) && obs_crs$epsg == 4326) {
    3857
  } else {
    4326
  }

  expect_error(
    utop_object(
      sf_obs,
      sf::st_transform(sf_pred, mismatch_crs),
      formula = "obs ~ 1"
    ),
    "different projections"
  )

  expect_error(
    utop_object(sf_obs, sf::st_set_crs(sf_pred, NA), formula = "obs ~ 1"),
    "only one of observations and prediction_locations have projection"
  )
})
