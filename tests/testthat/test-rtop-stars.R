test_that("stars vector cubes work end-to-end", {
  fixtures <- utop_stars_fixtures(n_obs = 8, n_pred = 4, n_time = 3)

  rtop_obj <- utop_object(
    fixtures$observations,
    fixtures$prediction_locations,
    formula = "obs ~ 1",
    params = utop_params(
      r_resol = 4,
      rs_type = "regular",
      debug_level = -1,
      nugget = FALSE
    )
  )

  expect_true(S7::S7_inherits(rtop_obj, Utop))
  expect_equal(unname(dim(rtop_obj@observations))[1], 8)
  expect_equal(unname(dim(rtop_obj@observations))[2], 3)

  rtop_obj <- utop_fit_variogram(rtop_obj, iprint = -1)
  expect_true(S7::S7_inherits(rtop_obj@variogram_model, UtopVariogramModel))

  result <- utop_krige(rtop_obj)
  expect_s3_class(result@predictions, "stars")
  expect_equal(unname(dim(result@predictions))[1], 4)
  expect_equal(unname(dim(result@predictions))[2], 3)
  expect_false(anyNA(result@predictions[["var1.pred"]]))
})

test_that("utop_variogram dispatches correctly for stars", {
  fixtures <- utop_stars_fixtures(n_obs = 6, n_pred = 3, n_time = 2)
  vario <- utop_variogram(fixtures$observations, formula = "obs ~ 1")

  expect_true(S7::S7_inherits(vario, UtopVariogram))
  expect_true("dist" %in% names(vario@data))
  expect_true("gamma" %in% names(vario@data))
})

test_that("stars pipeline tolerates nugget=TRUE", {
  fixtures <- utop_stars_fixtures(n_obs = 6, n_pred = 3, n_time = 2)

  rtop_obj <- utop_object(
    fixtures$observations,
    fixtures$prediction_locations,
    formula = "obs ~ 1",
    params = utop_params(
      r_resol = 4,
      rs_type = "regular",
      debug_level = -1,
      nugget = TRUE
    )
  )

  expect_true(S7::S7_inherits(rtop_obj, Utop))
  expect_true(!is.null(rtop_obj@overlap_obs))
  expect_error(utop_fit_variogram(rtop_obj, iprint = -1), NA)
})

test_that("stars time subsetting keeps the temporal dimension", {
  fixtures <- utop_stars_fixtures(n_obs = 5, n_pred = 2, n_time = 21)

  set.seed(1)
  expect_error(
    utop_object(
      fixtures$observations,
      fixtures$prediction_locations,
      formula = "obs ~ 1",
      params = utop_params(r_resol = 4, rs_type = "regular", debug_level = -1)
    ),
    NA
  )
  expect_error(
    utop_object(
      aperm(fixtures$observations, c("time", "geometry")),
      aperm(fixtures$prediction_locations, c("time", "geometry")),
      formula = "obs ~ 1",
      params = utop_params(r_resol = 4, rs_type = "regular", debug_level = -1)
    ),
    NA
  )
})

test_that("stars kriging supports lagged non-CV and CV paths", {
  fixtures <- utop_stars_fixtures(n_obs = 5, n_pred = 2, n_time = 4)
  rtop_obj <- utop_object(
    fixtures$observations,
    fixtures$prediction_locations,
    formula = "obs ~ 1",
    params = utop_params(
      r_resol = 4,
      rs_type = "regular",
      debug_level = -1,
      n_max = 3
    )
  )
  rtop_obj@variogram_model <- utop_variogram_model()

  pred <- utop_krige(
    rtop_obj,
    olags = rep(1, 5),
    plags = rep(0, 2),
    lag_exact = FALSE
  )
  expect_error(
    utop_krige(
      rtop_obj,
      cv = TRUE,
      olags = c(0, 0.25, 0.5, 0.75, 1),
      lag_exact = TRUE
    ),
    "subscript out of bounds|number of items to replace"
  )

  expect_s3_class(pred@predictions, "stars")
  expect_equal(unname(dim(pred@predictions)), c(2, 4))
  expect_true(!all(is.na(pred@predictions[["var1.pred"]])))
})
