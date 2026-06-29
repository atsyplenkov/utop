if (!identical(Sys.getenv("UTOP_DIAGNOSTICS"), "true")) {
  test_that("kriging and simulation diagnostics are enabled", {
    skip("Set UTOP_DIAGNOSTICS=true to run diagnostic integration tests")
  })
} else {
  spatial <- utop_spatial_subset_fixtures(n_obs = 12, n_pred = 2)
  unfit_base <- createRtopObject(
    spatial$observations,
    spatial$prediction_locations,
    formulaString = "obs ~ 1",
    params = modifyList(
      spatial$params,
      list(nugget = FALSE, model = "Ex1", gDist = FALSE)
    )
  )
  fit_base <- rtopFitVariogram(unfit_base, iprint = -1)

  spatial_small <- utop_spatial_subset_fixtures(n_obs = 8, n_pred = 2)
  fit_cloud <- rtopFitVariogram(
    createRtopObject(
      spatial_small$observations,
      spatial_small$prediction_locations,
      formulaString = "obs ~ 1",
      params = modifyList(
        spatial_small$params,
        list(nugget = FALSE, model = "Ex1", gDist = FALSE)
      )
    ),
    iprint = -1
  )

  test_that("rtopKrige covers sel, cvInfo, and uncertainty branches", {
    fit_sel <- fit_base
    fit_sel$predictionLocations <- fit_sel$predictionLocations[1, ]
    fit_sel$observations$unc <- rep(0.05, nrow(fit_sel$observations))
    fit_base$observations$unc <- rep(0.05, nrow(fit_base$observations))
    krige_sel <- rtopKrige(fit_sel, sel = 1)
    krige_cv <- rtopKrige(fit_base, cv = TRUE)

    expect_s3_class(krige_sel$predictions, "sf")
    expect_equal(nrow(krige_sel$predictions), 1)
    expect_true(all(
      c("var1.pred", "var1.var") %in% names(krige_sel$predictions)
    ))
    expect_true("cvInfo" %in% names(krige_cv))
    expect_true(nrow(krige_cv$cvInfo) > 0)
  })

  test_that("rtopSim covers error paths and the missing-area augmentation branch", {
    unfit <- createRtopObject(
      spatial$observations,
      spatial$prediction_locations,
      formulaString = "obs ~ 1",
      params = modifyList(spatial$params, list(nugget = FALSE, model = "Ex1"))
    )

    expect_error(
      rtopSim(unfit, nsim = 1, debug.level = -1),
      "Cannot do simulations without a variogram model"
    )

    fit <- rtopFitVariogram(unfit, iprint = -1)

    expect_error(
      rtopSim(fit, nsim = 1, replace = TRUE, debug.level = -1),
      "replaceNumber"
    )

    bad_replace <- fit
    bad_replace$predictionLocations <- bad_replace$observations[1:2, ]
    bad_replace$predictionLocations$replaceNumber <- c(
      1,
      nrow(bad_replace$observations) + 1
    )
    expect_error(
      rtopSim(bad_replace, nsim = 1, replace = TRUE, debug.level = -1),
      "does not correspond"
    )

    sim_input <- fit
    sim_input$predictionLocations <- sim_input$observations
    sim_input$predictionLocations$replaceNumber <- seq_len(nrow(
      sim_input$predictionLocations
    ))
    sim_input$predictionLocations$area <- NULL

    set.seed(1501)
    sim <- suppressWarnings(rtopSim(
      sim_input,
      nsim = 1,
      replace = FALSE,
      debug.level = -1
    ))

    expect_s3_class(sim$simulations, "sf")
    expect_true("area" %in% names(sim$simulations))
    expect_true("sim1" %in% names(sim$simulations))
  })

  test_that("checkVario and rtopCluster cover direct dispatch paths", {
    pdf_file <- tempfile(fileext = ".pdf")
    grDevices::pdf(pdf_file)
    on.exit(grDevices::dev.off(), add = TRUE)

    checked_rtop <- checkVario(
      fit_cloud,
      cloud = FALSE,
      gDist = FALSE,
      params = list(amul = 3, dmul = 3)
    )
    checked_rtop_cloud <- checkVario(
      fit_cloud,
      cloud = TRUE,
      gDist = FALSE,
      params = list(amul = 3, dmul = 3)
    )
    checked_model <- checkVario(
      fit_cloud$variogramModel,
      observations = fit_cloud$observations,
      sampleVariogram = fit_cloud$variogram,
      params = list(amul = 3, dmul = 3)
    )
    cloud_sample <- rtopVariogram(
      fit_cloud,
      params = modifyList(fit_cloud$params, list(cloud = TRUE, nugget = FALSE))
    )$variogramCloud
    checked_model_cloud <- checkVario(
      fit_cloud$variogramModel,
      observations = fit_cloud$observations,
      sampleVariogram = cloud_sample,
      params = list(amul = 3, dmul = 3)
    )

    cl <- rtopCluster(1, type = "PSOCK")
    on.exit(rtopCluster(1, action = "stop"), add = TRUE)
    cl_restart <- rtopCluster(1, action = "restart", type = "PSOCK")
    rtopCluster(1, action = "stop")

    expect_s3_class(checked_rtop, "rtop")
    expect_true("checkVario" %in% names(checked_rtop))
    expect_s3_class(checked_rtop_cloud, "rtop")
    expect_true("checkVario" %in% names(checked_rtop_cloud))
    expect_true(is.list(checked_model))
    expect_true(is.matrix(checked_model$vmats))
    expect_true(is.list(checked_model_cloud))
    expect_true(is.matrix(checked_model_cloud$vmats))
    expect_true(!is.null(cl))
    expect_true(!is.null(cl_restart))
  })
}
