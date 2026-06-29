fixtures <- utop_sf_subset_fixtures(n_obs = 12, n_pred = 5)

set.seed(1)
rtop_fitted <- createRtopObject(
  fixtures$observations,
  fixtures$prediction_locations,
  params = fixtures$params,
  formulaString = "obs ~1"
)
rtop_fitted <- rtopFitVariogram(rtop_fitted, iprint = -1)

test_that("sf kriging returns complete prediction fields", {
  rtop_pred <- rtopKrige(rtop_fitted)

  expect_s3_class(rtop_pred$predictions, "sf")
  expect_equal(nrow(rtop_pred$predictions), 5)
  expect_false(anyNA(rtop_pred$predictions$var1.pred))
  expect_false(anyNA(rtop_pred$predictions$var1.var))
})

test_that("sf kriging preserves the seeded prediction anchors", {
  rtop_pred <- rtopKrige(rtop_fitted)

  expect_lt(abs(rtop_pred$predictions$var1.pred[1] - 0.0107653735987488), 1e-7)
  expect_lt(abs(rtop_pred$predictions$var1.pred[2] - 0.0107530520520757), 1e-7)
})

test_that("sf cross-validation keeps the legacy correlation anchor", {
  rtop_cv <- rtopKrige(rtop_fitted, cv = TRUE)

  expect_equal(
    cor(rtop_cv$predictions$observed, rtop_cv$predictions$var1.pred),
    0.468766972817434,
    tolerance = 1e-7
  )
})
