test_that("stars vector cubes work end-to-end", {
  fixtures <- utop_stars_fixtures(n_obs = 8, n_pred = 4, n_time = 3)

  rtop_obj <- createRtopObject(
    fixtures$observations,
    fixtures$prediction_locations,
    formulaString = "obs ~ 1",
    params = list(
      rresol = 4,
      rstype = "regular",
      debug.level = -1,
      nugget = FALSE
    )
  )

  expect_s3_class(rtop_obj, "rtop")
  expect_equal(unname(dim(rtop_obj$observations))[1], 8)
  expect_equal(unname(dim(rtop_obj$observations))[2], 3)

  rtop_obj <- rtopFitVariogram(rtop_obj, iprint = -1)
  expect_s3_class(rtop_obj$variogramModel, "rtopVariogramModel")

  result <- rtopKrige(rtop_obj)
  expect_s3_class(result$predictions, "stars")
  expect_equal(unname(dim(result$predictions))[1], 4)
  expect_equal(unname(dim(result$predictions))[2], 3)
  expect_false(anyNA(result$predictions[["var1.pred"]]))
})

test_that("rtopVariogram dispatches correctly for stars", {
  fixtures <- utop_stars_fixtures(n_obs = 6, n_pred = 3, n_time = 2)
  vario <- rtopVariogram(fixtures$observations, formulaString = "obs ~ 1")

  expect_s3_class(vario, "rtopVariogram")
  expect_true("dist" %in% names(vario))
  expect_true("gamma" %in% names(vario))
})

test_that("stars pipeline tolerates nugget=TRUE", {
  fixtures <- utop_stars_fixtures(n_obs = 6, n_pred = 3, n_time = 2)

  rtop_obj <- createRtopObject(
    fixtures$observations,
    fixtures$prediction_locations,
    formulaString = "obs ~ 1",
    params = list(
      rresol = 4,
      rstype = "regular",
      debug.level = -1,
      nugget = TRUE
    )
  )

  expect_s3_class(rtop_obj, "rtop")
  expect_true("overlapObs" %in% names(rtop_obj))
  expect_error(rtopFitVariogram(rtop_obj, iprint = -1), NA)
})

test_that("stars time subsetting keeps the temporal dimension", {
  fixtures <- utop_stars_fixtures(n_obs = 5, n_pred = 2, n_time = 21)

  set.seed(1)
  expect_error(
    createRtopObject(
      fixtures$observations,
      fixtures$prediction_locations,
      formulaString = "obs ~ 1",
      params = list(rresol = 4, rstype = "regular", debug.level = -1)
    ),
    NA
  )
  expect_error(
    createRtopObject(
      aperm(fixtures$observations, c("time", "geometry")),
      aperm(fixtures$prediction_locations, c("time", "geometry")),
      formulaString = "obs ~ 1",
      params = list(rresol = 4, rstype = "regular", debug.level = -1)
    ),
    NA
  )
})

test_that("stars kriging supports lagged non-CV and CV paths", {
  fixtures <- utop_stars_fixtures(n_obs = 5, n_pred = 2, n_time = 4)
  rtop_obj <- createRtopObject(
    fixtures$observations,
    fixtures$prediction_locations,
    formulaString = "obs ~ 1",
    params = list(
      rresol = 4,
      rstype = "regular",
      debug.level = -1,
      nmax = 3
    )
  )
  rtop_obj$variogramModel <- rtopVariogramModel()

  pred <- rtopKrige(
    rtop_obj,
    olags = rep(1, 5),
    plags = rep(0, 2),
    lagExact = FALSE
  )
  cv <- rtopKrige(
    rtop_obj,
    cv = TRUE,
    olags = c(0, 0.25, 0.5, 0.75, 1),
    lagExact = TRUE
  )

  expect_s3_class(pred$predictions, "stars")
  expect_s3_class(cv$predictions, "stars")
  expect_equal(unname(dim(pred$predictions)), c(2, 4))
  expect_equal(unname(dim(cv$predictions)), c(5, 4))
  expect_true(any(!is.na(pred$predictions[["var1.pred"]])))
  expect_true(any(!is.na(cv$predictions[["var1.pred"]])))
})
