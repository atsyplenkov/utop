spatial <- utop_spatial_subset_fixtures(n_obs = 8, n_pred = 2)
sf_fixtures <- utop_sf_subset_fixtures(n_obs = 6, n_pred = 2)

sp_obs <- spatial$observations
sp_pred <- spatial$prediction_locations
sf_obs <- sf_fixtures$observations

params_sp <- modifyList(spatial$params, list(nugget = FALSE, model = "Ex1"))
params_sf_cloud <- modifyList(
  sf_fixtures$params,
  list(cloud = TRUE, nugget = FALSE, model = "Ex1")
)
params_sf_cloud_false <- modifyList(
  sf_fixtures$params,
  list(
    cloud = FALSE,
    nugget = FALSE,
    model = "Ex1",
    g_dist_est = FALSE,
    g_dist_pred = FALSE
  )
)

vario_sp <- utop_variogram(sp_obs, formula = "obs ~ 1", params = params_sp)
vario_sp_cloud <- utop_variogram(
  sp_obs,
  formula = "obs ~ 1",
  params = modifyList(params_sp, list(cloud = TRUE))
)
vario_sf <- utop_variogram(
  sf_obs,
  formula = "obs ~ 1",
  params = params_sf_cloud_false
)
vario_sf_cloud <- utop_variogram(
  sf_obs,
  formula = "obs ~ 1",
  params = params_sf_cloud
)

disc_sp <- utop_disc(
  sp_obs,
  params = list(rs_type = "regular", r_resol = 4, debug_level = -1)
)
disc_sf <- utop_disc(
  sf_obs,
  params = list(rs_type = "regular", r_resol = 4, debug_level = -1)
)

fit_sp <- utop_fit_variogram(
  utop_object(sp_obs, sp_pred, formula = "obs ~ 1", params = params_sp),
  iprint = -1
)
fit_sf <- utop_fit_variogram(
  utop_object(
    sf_obs,
    sf_fixtures$prediction_locations,
    formula = "obs ~ 1",
    params = params_sf_cloud
  ),
  iprint = -1
)

varmat_cached <- utop_var_mat(fit_sp)
vm_default <- utop_variogram_model()

test_that("variogram model construction and updates cover default and mutation branches", {
  vm_default_local <- utop_variogram_model()
  expect_true(S7::S7_inherits(vm_default_local, UtopVariogramModel))
  expect_identical(vm_default_local@model, "Ex1")
  expect_length(vm_default_local@params, 5)

  vm_from_obs <- utop_variogram_model(
    observations = sp_obs,
    formula = as.formula("obs ~ 1")
  )
  expect_true(S7::S7_inherits(vm_from_obs, UtopVariogramModel))
  expect_length(vm_from_obs@params, 5)

  vm_mult <- utop_update_variogram(
    vm_default_local,
    sill = 2,
    range = 3,
    nugget = 4,
    exp = 5,
    exp0 = 6,
    action = "mult"
  )
  vm_replace <- utop_update_variogram(
    vm_default_local,
    sill = 99,
    action = "replace"
  )
  vm_add <- utop_update_variogram(vm_default_local, exp = 1.5, action = "add")

  expect_equal(vm_mult@params[1], vm_default_local@params[1] * 2)
  expect_equal(vm_mult@params[2], vm_default_local@params[2] * 3)
  expect_equal(vm_replace@params[1], 99)
  expect_equal(vm_add@params[4], vm_default_local@params[4] * 1.5)
})

test_that("variogram fitting, discretization, and matrices cover spatial and sf paths", {
  expect_true(S7::S7_inherits(vario_sp, UtopVariogram))
  expect_true(S7::S7_inherits(vario_sp_cloud, UtopVariogramCloud))
  expect_true(S7::S7_inherits(vario_sf, UtopVariogram))
  expect_true(S7::S7_inherits(vario_sf_cloud, UtopVariogramCloud))

  disc_from_vario <- utop_disc(
    vario_sp,
    params = list(h_resol = 2, hs_type = "regular", debug_level = -1)
  )
  expect_length(disc_from_vario, nrow(vario_sp@data))
  expect_length(disc_from_vario[[1]], 2)

  expect_length(disc_sp, nrow(sp_obs))
  expect_length(disc_sf, nrow(sf_obs))

  gd <- utop_g_dist(disc_sp, params = list(debug_level = -1))
  expect_equal(dim(gd), c(length(disc_sp), length(disc_sp)))
  expect_equal(gd, t(gd))
  expect_true(all(is.finite(diag(gd))))

  expect_true(S7::S7_inherits(fit_sp@variogram_model, UtopVariogramModel))
  expect_true(S7::S7_inherits(fit_sf@variogram_model, UtopVariogramModel))
  expect_length(fit_sp@variogram_model@params, 5)
  expect_length(fit_sf@variogram_model@params, 5)
  expect_identical(fit_sp@variogram_model@params[3], 0)
  expect_true(!is.null(fit_sf@variogram_cloud))

  varmat_reuse <- utop_var_mat(varmat_cached)
  varmat_cv <- utop_var_mat(varmat_cached, params = list(cv = TRUE))

  varmat_list <- utop_var_mat(
    sp_obs,
    sp_pred,
    variogram_model = vm_default,
    g_dist_pred = FALSE,
    nugget = FALSE,
    debug_level = -1
  )

  vmod_nugget <- utop_variogram_model(nugget = 0.1)
  varmat_no_nugget <- utop_var_mat(
    sp_obs,
    sp_pred,
    variogram_model = vmod_nugget,
    nugget = FALSE,
    r_resol = 4,
    rs_type = "regular",
    debug_level = -1
  )
  expect_error(
    varmat_nugget <- utop_var_mat(
      sp_obs,
      sp_pred,
      variogram_model = vmod_nugget,
      nugget = TRUE,
      r_resol = 4,
      rs_type = "regular",
      debug_level = -1
    ),
    NA
  )
  expect_error(
    varmat_obs_nugget <- utop_var_mat(
      sp_obs,
      variogram_model = vmod_nugget,
      nugget = TRUE,
      r_resol = 4,
      rs_type = "regular",
      debug_level = -1
    ),
    NA
  )

  expect_equal(dim(varmat_cached@var_mat_obs), c(nrow(sp_obs), nrow(sp_obs)))
  expect_equal(varmat_cached@var_mat_obs, t(varmat_cached@var_mat_obs))
  expect_equal(
    dim(varmat_cached@var_mat_pred_obs),
    c(nrow(sp_obs), nrow(sp_pred))
  )
  expect_equal(dim(varmat_cached@var_mat_pred), c(nrow(sp_pred), 1))
  expect_true(all(is.finite(diag(varmat_cached@var_mat_obs))))

  expect_equal(dim(varmat_list@var_mat_obs), c(nrow(sp_obs), nrow(sp_obs)))
  expect_equal(
    dim(varmat_list@var_mat_pred_obs),
    c(nrow(sp_obs), nrow(sp_pred))
  )
  expect_true(!is.null(varmat_list@var_mat_pred_obs))
  expect_equal(dim(varmat_nugget@var_mat_obs), c(nrow(sp_obs), nrow(sp_obs)))
  expect_equal(
    dim(varmat_nugget@var_mat_pred_obs),
    c(nrow(sp_obs), nrow(sp_pred))
  )
  expect_equal(dim(varmat_obs_nugget), c(nrow(sp_obs), nrow(sp_obs)))
  expect_true(all(is.finite(varmat_obs_nugget)))
  expect_true(any(
    varmat_nugget@var_mat_pred_obs > varmat_no_nugget@var_mat_pred_obs
  ))

  expect_identical(varmat_cached@var_mat_obs, varmat_reuse@var_mat_obs)
  expect_identical(
    varmat_cached@var_mat_pred_obs,
    varmat_reuse@var_mat_pred_obs
  )
  expect_identical(varmat_cached@var_mat_obs, varmat_cv@var_mat_obs)

  expect_true(all(is.finite(varmat_cached@var_mat_obs)))
  expect_true(all(is.finite(varmat_list@var_mat_pred_obs)))
})

test_that("overlap helpers cover their direct branches", {
  overlap_self <- findOverlap(sp_obs[1:4, ], debug.level = 0)
  overlap_cross <- findOverlap(sp_obs[1:4, ], sp_obs[5:6, ], debug.level = 0)

  expect_equal(dim(overlap_self), c(4, 4))
  expect_equal(overlap_self, t(overlap_self))
  expect_true(all(diag(overlap_self) > 0))
  expect_equal(dim(overlap_cross), c(4, 2))
  expect_true(all(overlap_cross >= 0))
})
