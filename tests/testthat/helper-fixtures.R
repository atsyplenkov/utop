utop_extdata_path <- function() {
  extdata_path <- system.file("extdata", package = "utop")
  if (nzchar(extdata_path)) {
    return(extdata_path)
  }

  if (requireNamespace("testthat", quietly = TRUE)) {
    extdata_path <- testthat::test_path("..", "..", "inst", "extdata")
    if (dir.exists(extdata_path)) {
      return(normalizePath(extdata_path))
    }
  }

  stop("Could not locate utop extdata for tests")
}

utop_read_sf <- function(layer) {
  sf::st_read(utop_extdata_path(), layer, quiet = TRUE)
}

utop_demo_data <- function() {
  readRDS(file.path(utop_extdata_path(), "demo.rds"))
}

utop_test_params <- function(...) {
  defaults <- list(
    g_dist_est = TRUE,
    g_dist_pred = TRUE,
    cloud = FALSE,
    r_resol = 25,
    h_resol = 3,
    debug_level = -1
  )
  do.call(utop_params, utils::modifyList(defaults, list(...)))
}

utop_update_params <- function(params, ...) {
  update_utop_params(params, list(...))
}


utop_spatial_fixtures <- function() {
  observations <- utop_read_sf("observations")
  prediction_locations <- utop_read_sf("predictionLocations")

  observations <- observations[1:30, ]
  prediction_locations <- prediction_locations[1:2, ]
  observations$obs <- observations$QSUMMER_OB / observations$AREASQKM

  list(
    observations = observations,
    prediction_locations = prediction_locations,
    params = utop_test_params()
  )
}


utop_sf_subset_fixtures <- function(
  n_obs = 10,
  n_pred = 5,
  params = utop_test_params()
) {
  observations <- utop_read_sf("observations")
  prediction_locations <- utop_read_sf("predictionLocations")
  observations$obs <- observations$QSUMMER_OB / observations$AREASQKM

  list(
    observations = observations[seq_len(n_obs), ],
    prediction_locations = prediction_locations[seq_len(n_pred), ],
    params = params
  )
}

utop_spatial_subset_fixtures <- function(
  n_obs = 10,
  n_pred = 2,
  params = utop_test_params()
) {
  fixtures <- utop_spatial_fixtures()

  list(
    observations = fixtures$observations[seq_len(n_obs), ],
    prediction_locations = fixtures$prediction_locations[seq_len(n_pred), ],
    params = params
  )
}
