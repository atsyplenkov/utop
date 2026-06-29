fixtures <- utop_spatial_fixtures()

set.seed(1501)
rtop_fitted <- createRtopObject(
  fixtures$observations,
  fixtures$prediction_locations,
  params = fixtures$params,
  formulaString = "obs ~ 1"
)
rtop_fitted <- rtopFitVariogram(rtop_fitted, iprint = -1)

test_that("kriging returns expected prediction structures", {
  rtop_cv <- rtopKrige(rtop_fitted, cv = TRUE)
  rtop_pred <- rtopKrige(rtop_fitted)

  expect_s3_class(rtop_fitted$variogramModel, "rtopVariogramModel")

  expect_s3_class(rtop_cv$predictions, "sf")
  expect_s3_class(rtop_pred$predictions, "sf")
  expect_equal(nrow(rtop_cv$predictions), 30)
  expect_equal(nrow(rtop_pred$predictions), 2)
  expect_true(all(
    c("observed", "var1.pred", "var1.var") %in% names(rtop_cv$predictions)
  ))
  expect_true(all(c("var1.pred", "var1.var") %in% names(rtop_pred$predictions)))
})

test_that("kriging reuses the legacy semivariance path consistently", {
  rtop_cv <- rtopKrige(rtop_fitted, cv = TRUE)
  rtop_pred <- rtopKrige(rtop_fitted)
  varmat <- varMat(
    fixtures$observations,
    fixtures$prediction_locations,
    variogramModel = rtop_fitted$variogramModel,
    gDistEst = TRUE,
    gDistPred = TRUE,
    rresol = 25,
    hresol = 3
  )
  rtop_reuse <- rtopKrige(rtop_cv)

  expect_s3_class(rtop_reuse$predictions, "sf")
  expect_equal(nrow(rtop_reuse$predictions), 2)
  expect_true(isTRUE(all.equal(varmat$varMatObs, rtop_cv$varMatObs)))
  expect_true(isTRUE(all.equal(rtop_reuse$predictions, rtop_pred$predictions)))
})

test_that("spatial cross-validation keeps the legacy correlation anchor", {
  rtop_cv <- rtopKrige(rtop_fitted, cv = TRUE)

  expect_equal(
    cor(rtop_cv$predictions$observed, rtop_cv$predictions$var1.pred),
    0.165756902165012,
    tolerance = 1e-7
  )
})

test_that("spatial variogram updates rebuild semivariance matrices", {
  rtop_reuse <- rtopKrige(rtopKrige(rtop_fitted, cv = TRUE))

  rtop_updated <- varMat(rtop_reuse)
  rtop_updated <- updateRtopVariogram(rtop_updated, exp = 1.5, action = "mult")
  rtop_updated_mat <- varMat(rtop_updated)

  expect_false(isTRUE(all.equal(
    rtop_updated$varMatObs,
    rtop_updated_mat$varMatObs
  )))
  expect_true(!is.null(rtop_updated_mat$varMatObs))
})

test_that("spatial simulation stays anchored to the seeded legacy run", {
  set.seed(1501)
  rtop_sim_2 <- rtopSim(rtop_fitted, nsim = 2, logdist = TRUE, debug.level = -1)
  rtop_sim_input <- rtop_fitted
  rtop_sim_input$predictionLocations <- rtop_sim_input$observations
  rtop_sim_input$observations$unc <- var(rtop_sim_2$observations$obs) *
    min(rtop_sim_2$observations$area) /
    rtop_sim_2$observations$area
  rtop_sim_input$predictionLocations$replaceNumber <- seq_len(nrow(
    rtop_sim_input$predictionLocations
  ))
  rtop_sim_3 <- rtopSim(
    rtop_sim_input,
    nsim = 3,
    replace = TRUE,
    debug.level = -1
  )

  expect_equal(
    rtop_sim_2$simulations$sim1[1],
    0.0118082328747604,
    tolerance = 1e-7
  )
  expect_equal(
    rtop_sim_2$simulations$sim2[1],
    0.0114985387672817,
    tolerance = 1e-7
  )
  expect_equal(
    rtop_sim_2$simulations$sim2[2],
    0.0103816065249376,
    tolerance = 1e-7
  )
  expect_equal(
    rtop_sim_3$simulations$sim1[1],
    0.0126686566529176,
    tolerance = 1e-7
  )
  expect_equal(
    rtop_sim_3$simulations$sim2[1],
    0.0120647228006132,
    tolerance = 1e-7
  )
  expect_equal(
    rtop_sim_3$simulations$sim3[14],
    0.019830061178791,
    tolerance = 1e-7
  )
})

