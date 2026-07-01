fixtures <- utop_sf_subset_fixtures(n_obs = 12, n_pred = 5)

set.seed(1)
rtop_fitted <- utop_object(
  fixtures$observations,
  fixtures$prediction_locations,
  params = fixtures$params,
  formula = "obs ~1"
)
rtop_fitted <- utop_fit_variogram(rtop_fitted, iprint = -1)

test_that("sf kriging returns complete prediction fields", {
  rtop_pred <- utop_krige(rtop_fitted)

  expect_s3_class(rtop_pred@predictions, "sf")
  expect_equal(nrow(rtop_pred@predictions), 5)
  expect_false(anyNA(rtop_pred@predictions$var1.pred))
  expect_false(anyNA(rtop_pred@predictions$var1.var))
})
