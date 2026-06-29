# jarl-ignore-file internal_function: testing internal functions

spatial <- utop_spatial_subset_fixtures(n_obs = 8, n_pred = 2)
sf_fixtures <- utop_sf_subset_fixtures(n_obs = 6, n_pred = 2)

sp_obs <- spatial$observations
sp_pred <- spatial$prediction_locations
sf_obs <- sf_fixtures$observations

params_sp <- modifyList(
  spatial$params,
  list(nugget = FALSE, model = "Ex1", gDist = TRUE)
)
params_sf_cloud <- modifyList(
  sf_fixtures$params,
  list(cloud = TRUE, nugget = FALSE, model = "Ex1", gDist = TRUE)
)
params_sf_cloud_false <- modifyList(
  sf_fixtures$params,
  list(cloud = FALSE, nugget = FALSE, model = "Ex1", gDist = FALSE)
)

vario_sp <- rtopVariogram(sp_obs, formulaString = "obs ~ 1", params = params_sp)
vario_sp_cloud <- rtopVariogram(
  sp_obs,
  formulaString = "obs ~ 1",
  params = modifyList(params_sp, list(cloud = TRUE))
)
vario_sf <- rtopVariogram(
  sf_obs,
  formulaString = "obs ~ 1",
  params = params_sf_cloud_false
)
vario_sf_cloud <- rtopVariogram(
  sf_obs,
  formulaString = "obs ~ 1",
  params = params_sf_cloud
)

disc_sp <- rtopDisc(
  sp_obs,
  params = list(rstype = "regular", rresol = 4, debug.level = -1)
)
disc_sf <- rtopDisc(
  sf_obs,
  params = list(rstype = "regular", rresol = 4, debug.level = -1)
)

fit_sp <- rtopFitVariogram(
  createRtopObject(
    sp_obs,
    sp_pred,
    formulaString = "obs ~ 1",
    params = params_sp
  ),
  iprint = -1
)
fit_sf <- rtopFitVariogram(
  createRtopObject(
    sf_obs,
    sf_fixtures$prediction_locations,
    formulaString = "obs ~ 1",
    params = params_sf_cloud
  ),
  iprint = -1
)

varmat_cached <- varMat(fit_sp)
vm_default <- rtopVariogramModel()

test_that("variogram model construction and updates cover default and mutation branches", {
  vm_default_local <- rtopVariogramModel()
  expect_s3_class(vm_default_local, "rtopVariogramModel")
  expect_identical(vm_default_local$model, "Ex1")
  expect_length(vm_default_local$params, 5)

  vm_from_obs <- rtopVariogramModel(
    observations = sp_obs,
    formulaString = as.formula("obs ~ 1")
  )
  expect_s3_class(vm_from_obs, "rtopVariogramModel")
  expect_length(vm_from_obs$params, 5)

  vm_mult <- updateRtopVariogram(
    vm_default_local,
    sill = 2,
    range = 3,
    nugget = 4,
    exp = 5,
    exp0 = 6,
    action = "mult"
  )
  vm_replace <- updateRtopVariogram(
    vm_default_local,
    sill = 99,
    action = "replace"
  )
  vm_add <- updateRtopVariogram(vm_default_local, exp = 1.5, action = "add")

  expect_equal(vm_mult$params[1], vm_default_local$params[1] * 2)
  expect_equal(vm_mult$params[2], vm_default_local$params[2] * 3)
  expect_equal(vm_replace$params[1], 99)
  expect_equal(vm_add$params[4], vm_default_local$params[4] * 1.5)
})

test_that("variogram fitting, discretization, and matrices cover spatial and sf paths", {
  expect_s3_class(vario_sp, "rtopVariogram")
  expect_s3_class(vario_sp_cloud, "rtopVariogramCloud")
  expect_s3_class(vario_sf, "rtopVariogram")
  expect_s3_class(vario_sf_cloud, "rtopVariogramCloud")

  disc_from_vario <- rtopDisc(
    vario_sp,
    params = list(hresol = 2, hstype = "regular", debug.level = -1)
  )
  expect_length(disc_from_vario, nrow(vario_sp))
  expect_length(disc_from_vario[[1]], 2)

  expect_length(disc_sp, nrow(sp_obs))
  expect_length(disc_sf, nrow(sf_obs))

  gd <- gDist(disc_sp, params = list(debug.level = -1))
  expect_equal(dim(gd), c(length(disc_sp), length(disc_sp)))
  expect_equal(gd, t(gd))
  expect_true(all(is.finite(diag(gd))))

  expect_s3_class(fit_sp$variogramModel, "rtopVariogramModel")
  expect_s3_class(fit_sf$variogramModel, "rtopVariogramModel")
  expect_length(fit_sp$variogramModel$params, 5)
  expect_length(fit_sf$variogramModel$params, 5)
  expect_identical(fit_sp$variogramModel$params[3], 0)
  expect_true(!is.null(fit_sf$variogramCloud))

  varmat_reuse <- varMat(varmat_cached)
  varmat_cv <- varMat(varmat_cached, params = list(cv = TRUE))

  varmat_list <- varMat(
    sp_obs,
    sp_pred,
    variogramModel = vm_default,
    gDistPred = FALSE,
    nugget = FALSE,
    debug.level = -1
  )

  vmod_nugget <- rtopVariogramModel(nugget = 0.1)
  varmat_no_nugget <- varMat(
    sp_obs,
    sp_pred,
    variogramModel = vmod_nugget,
    nugget = FALSE,
    rresol = 4,
    rstype = "regular",
    debug.level = -1
  )
  expect_error(
    varmat_nugget <- varMat(
      sp_obs,
      sp_pred,
      variogramModel = vmod_nugget,
      nugget = TRUE,
      rresol = 4,
      rstype = "regular",
      debug.level = -1
    ),
    NA
  )
  expect_error(
    varmat_obs_nugget <- varMat(
      sp_obs,
      variogramModel = vmod_nugget,
      nugget = TRUE,
      rresol = 4,
      rstype = "regular",
      debug.level = -1
    ),
    NA
  )

  expect_equal(dim(varmat_cached$varMatObs), c(nrow(sp_obs), nrow(sp_obs)))
  expect_equal(varmat_cached$varMatObs, t(varmat_cached$varMatObs))
  expect_equal(dim(varmat_cached$varMatPredObs), c(nrow(sp_obs), nrow(sp_pred)))
  expect_equal(dim(varmat_cached$varMatPred), c(nrow(sp_pred), 1))
  expect_true(all(is.finite(diag(varmat_cached$varMatObs))))

  expect_equal(dim(varmat_list$varMatObs), c(nrow(sp_obs), nrow(sp_obs)))
  expect_equal(dim(varmat_list$varMatPredObs), c(nrow(sp_obs), nrow(sp_pred)))
  expect_true(!is.null(varmat_list$varMatPredObs))
  expect_equal(dim(varmat_nugget$varMatObs), c(nrow(sp_obs), nrow(sp_obs)))
  expect_equal(dim(varmat_nugget$varMatPredObs), c(nrow(sp_obs), nrow(sp_pred)))
  expect_equal(dim(varmat_obs_nugget), c(nrow(sp_obs), nrow(sp_obs)))
  expect_true(all(is.finite(varmat_obs_nugget)))
  expect_true(any(
    varmat_nugget$varMatPredObs > varmat_no_nugget$varMatPredObs
  ))

  expect_identical(varmat_cached$varMatObs, varmat_reuse$varMatObs)
  expect_identical(varmat_cached$varMatPredObs, varmat_reuse$varMatPredObs)
  expect_identical(varmat_cached$varMatObs, varmat_cv$varMatObs)

  expect_true(all(is.finite(varmat_cached$varMatObs)))
  expect_true(all(is.finite(varmat_list$varMatPredObs)))
})

test_that("overlap and wrapper helpers cover their direct branches", {
  overlap_self <- findOverlap(sp_obs[1:4, ], debug.level = 0)
  overlap_cross <- findOverlap(sp_obs[1:4, ], sp_obs[5:6, ], debug.level = 0)

  expect_equal(dim(overlap_self), c(4, 4))
  expect_equal(overlap_self, t(overlap_self))
  expect_true(all(diag(overlap_self) > 0))
  expect_equal(dim(overlap_cross), c(4, 2))
  expect_true(all(overlap_cross >= 0))

  est <- utop:::estimateParameters.rtop(createRtopObject(
    sp_obs,
    spatial$prediction_locations,
    formulaString = "obs ~ 1",
    params = modifyList(spatial$params, list(nugget = FALSE, gDist = FALSE))
  ))
  mp <- utop:::methodParameters.rtop(est)

  expect_s3_class(est, "rtop")
  expect_s3_class(est$variogramModel, "rtopVariogramModel")
  expect_match(mp$methodParameters, "vmodel")
})
