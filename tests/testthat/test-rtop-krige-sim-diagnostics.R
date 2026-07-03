if (!identical(Sys.getenv("UTOP_DIAGNOSTICS"), "true")) {
  test_that("kriging and simulation diagnostics are enabled", {
    skip("Set UTOP_DIAGNOSTICS=true to run diagnostic integration tests")
  })
} else {
  spatial <- utop_spatial_subset_fixtures(n_obs = 12, n_pred = 2)
  unfit_base <- utop_object(
    spatial$observations,
    spatial$prediction_locations,
    formula = "obs ~ 1",
    params = utop_update_params(
      spatial$params,
      nugget = FALSE,
      model = "Ex1",
      g_dist_est = FALSE,
      g_dist_pred = FALSE
    )
  )
  fit_base <- utop_fit_variogram(unfit_base, iprint = -1)

  spatial_small <- utop_spatial_subset_fixtures(n_obs = 8, n_pred = 2)
  fit_cloud <- utop_fit_variogram(
    utop_object(
      spatial_small$observations,
      spatial_small$prediction_locations,
      formula = "obs ~ 1",
      params = utop_update_params(
        spatial_small$params,
        nugget = FALSE,
        model = "Ex1",
        g_dist_est = FALSE,
        g_dist_pred = FALSE
      )
    ),
    iprint = -1
  )

  test_that("utop_krige covers sel, cvInfo, and uncertainty branches", {
    fit_sel <- fit_base
    fit_sel$prediction_locations <- fit_sel$prediction_locations[1, ]
    fit_sel@observations$unc <- rep(0.05, nrow(fit_sel@observations))
    fit_base@observations$unc <- rep(0.05, nrow(fit_base@observations))
    krige_sel <- utop_krige(fit_sel, sel = 1)
    krige_cv <- utop_krige(fit_base, cv = TRUE)

    expect_s3_class(krige_sel@predictions, "sf")
    expect_equal(nrow(krige_sel@predictions), 1)
    expect_true(all(
      c("var1.pred", "var1.var") %in% names(krige_sel@predictions)
    ))
    expect_true("cvInfo" %in% names(krige_cv))
    expect_true(nrow(krige_cv@cv_info) > 0)
  })

  test_that("utop_sim covers error paths and the missing-area augmentation branch", {
    unfit <- utop_object(
      spatial$observations,
      spatial$prediction_locations,
      formula = "obs ~ 1",
      params = utop_update_params(spatial$params, nugget = FALSE, model = "Ex1")
    )

    expect_error(
      utop_sim(unfit, nsim = 1, debug_level = -1),
      "Cannot do simulations without a variogram model"
    )

    fit <- utop_fit_variogram(unfit, iprint = -1)

    expect_error(
      utop_sim(fit, nsim = 1, replace = TRUE, debug_level = -1),
      "replaceNumber"
    )

    bad_replace <- fit
    bad_replace$prediction_locations <- bad_replace@observations[1:2, ]
    bad_replace$prediction_locations$replaceNumber <- c(
      1,
      nrow(bad_replace@observations) + 1
    )
    expect_error(
      utop_sim(bad_replace, nsim = 1, replace = TRUE, debug_level = -1),
      "does not correspond"
    )

    sim_input <- fit
    sim_input$prediction_locations <- sim_input@observations
    sim_input$prediction_locations$replaceNumber <- seq_len(nrow(
      sim_input$prediction_locations
    ))
    sim_input$prediction_locations$area <- NULL

    set.seed(1501)
    sim <- suppressWarnings(utop_sim(
      sim_input,
      nsim = 1,
      replace = FALSE,
      debug_level = -1
    ))

    expect_s3_class(sim@simulations, "sf")
    expect_true("area" %in% names(sim@simulations))
    expect_true("sim1" %in% names(sim@simulations))
  })

  test_that("utop_check_variogram and utop_cluster cover direct dispatch paths", {
    pdf_file <- tempfile(fileext = ".pdf")
    grDevices::pdf(pdf_file)
    on.exit(grDevices::dev.off(), add = TRUE)

    checked_rtop <- utop_check_variogram(
      fit_cloud,
      cloud = FALSE,
      gDist = FALSE,
      params = utop_params(amul = 3, dmul = 3)
    )
    checked_rtop_cloud <- utop_check_variogram(
      fit_cloud,
      cloud = TRUE,
      gDist = FALSE,
      params = utop_params(amul = 3, dmul = 3)
    )
    checked_model <- utop_check_variogram(
      fit_cloud@variogram_model,
      observations = fit_cloud@observations,
      sample_variogram = fit_cloud@variogram,
      params = utop_params(amul = 3, dmul = 3)
    )
    cloud_sample <- utop_variogram(
      fit_cloud,
      params = utop_update_params(
        fit_cloud@params,
        cloud = TRUE,
        nugget = FALSE
      )
    )@variogram_cloud
    checked_model_cloud <- utop_check_variogram(
      fit_cloud@variogram_model,
      observations = fit_cloud@observations,
      sample_variogram = cloud_sample,
      params = utop_params(amul = 3, dmul = 3)
    )

    cl <- utop_cluster(1, type = "PSOCK")
    on.exit(utop_cluster(1, action = "stop"), add = TRUE)
    cl_restart <- utop_cluster(1, action = "restart", type = "PSOCK")
    utop_cluster(1, action = "stop")

    expect_true(S7::S7_inherits(checked_rtop, Utop))
    expect_true("utop_check_variogram" %in% names(checked_rtop))
    expect_true(S7::S7_inherits(checked_rtop_cloud, Utop))
    expect_true("utop_check_variogram" %in% names(checked_rtop_cloud))
    expect_true(is.list(checked_model))
    expect_true(is.matrix(checked_model$vmats))
    expect_true(is.list(checked_model_cloud))
    expect_true(is.matrix(checked_model_cloud$vmats))
    expect_true(!is.null(cl))
    expect_true(!is.null(cl_restart))
  })
}
