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

  expect_s4_class(rtop_cv$predictions, "SpatialPolygonsDataFrame")
  expect_s4_class(rtop_pred$predictions, "SpatialPolygonsDataFrame")
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

  expect_s4_class(rtop_reuse$predictions, "SpatialPolygonsDataFrame")
  expect_equal(nrow(rtop_reuse$predictions), 2)
  expect_true(isTRUE(all.equal(varmat$varMatObs, rtop_cv$varMatObs)))
  expect_true(isTRUE(all.equal(rtop_reuse$predictions, rtop_pred$predictions)))
})

test_that("spatial cross-validation keeps the legacy correlation anchor", {
  rtop_cv <- rtopKrige(rtop_fitted, cv = TRUE)

  expect_equal(
    cor(rtop_cv$predictions$observed, rtop_cv$predictions$var1.pred),
    0.167324265652400,
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
    rtop_sim_2$simulations@data$sim1[1],
    0.0118349997842537,
    tolerance = 1e-7
  )
  expect_equal(
    rtop_sim_2$simulations@data$sim2[1],
    0.0115152470755103,
    tolerance = 1e-7
  )
  expect_equal(
    rtop_sim_2$simulations@data$sim2[2],
    0.0103075931751047,
    tolerance = 1e-7
  )
  expect_equal(
    rtop_sim_3$simulations@data$sim1[1],
    0.0125934486588040,
    tolerance = 1e-7
  )
  expect_equal(
    rtop_sim_3$simulations@data$sim2[1],
    0.0120003913823317,
    tolerance = 1e-7
  )
  expect_equal(
    rtop_sim_3$simulations@data$sim3[14],
    0.0199191982053638,
    tolerance = 1e-7
  )
})

test_that("intamap interpolation matches direct spatial kriging", {
  skip_if_not_installed("intamap")
  probe <- try(useRtopWithIntamap(), silent = TRUE)
  if (inherits(probe, "try-error")) {
    skip("intamap does not expose estimateParameters() on this system")
  }

  set.seed(1501)
  rtop_pred <- rtopKrige(rtop_fitted)

  output <- interpolate(
    fixtures$observations,
    fixtures$prediction_locations,
    optList = list(
      formulaString = obs ~ 1,
      gDist = TRUE,
      cloud = FALSE,
      nmax = 10,
      rresol = 25,
      hresol = 3
    ),
    methodName = "rtop",
    iprint = -1
  )

  expect_true(isTRUE(all.equal(
    rtop_pred$predictions@data$var1.pred,
    output$predictions@data$var1.pred
  )))
  expect_true(isTRUE(all.equal(
    rtop_pred$predictions@data$var1.var,
    output$predictions@data$var1.var
  )))
})
