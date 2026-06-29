utop_stars_cube <- function(space, time, values) {
  grid <- expand.grid(
    space = seq_len(nrow(space)),
    time = seq_along(time),
    KEEP.OUT.ATTRS = FALSE
  )
  long <- sf::st_drop_geometry(space)[grid$space, , drop = FALSE]
  long$time <- time[grid$time]
  for (nm in names(values)) {
    long[[nm]] <- as.vector(values[[nm]])
  }
  sf::st_geometry(long) <- sf::st_geometry(space)[grid$space]
  stars::st_as_stars(sf::st_as_sf(long), dims = c("geometry", "time"))
}

utop_stars_fixtures <- function(n_obs = 8, n_pred = 4, n_time = 3) {
  set.seed(42)
  x <- c(0, 2, 2, 0, 0)
  y <- c(0, 0, 2, 2, 0)

  make_polys <- function(n, prefix) {
    geoms <- lapply(seq_len(n), function(i) {
      px <- x + runif(1, 0, 20)
      py <- y + runif(1, 0, 20)
      sf::st_polygon(list(cbind(px, py)))
    })
    sfobj <- sf::st_sf(
      id = paste0(prefix, seq_len(n)),
      geometry = sf::st_sfc(geoms)
    )
    sfobj$area <- as.numeric(sf::st_area(sfobj))
    sfobj
  }

  sf_obs <- make_polys(n_obs, "P")
  sf_pred <- make_polys(n_pred, "Pr")
  time <- as.POSIXct("2020-01-01", tz = "UTC") + seq_len(n_time) * 86400

  base_obs <- runif(n_obs, 1, 10)
  obs <- matrix(NA_real_, nrow = n_obs, ncol = n_time)
  for (i in seq_len(n_obs)) {
    obs[i, ] <- base_obs[i] + rnorm(n_time, 0, 0.5)
  }
  pred <- matrix(0, nrow = n_pred, ncol = n_time)

  list(
    observations = utop_stars_cube(sf_obs, time, list(obs = obs)),
    prediction_locations = utop_stars_cube(sf_pred, time, list(var1 = pred)),
    time = time
  )
}

utop_tndtk_stars_fixture <- function(
  n_obs = 8,
  n_pred = 3,
  dates = as.Date("2006-01-01") + 0:4
) {
  demo <- utop_demo_data()
  observations <- demo$gauged_catchments[seq_len(n_obs), ]
  prediction_locations <- demo$ungauged_catchments[seq_len(n_pred), ]
  observations$area <- observations$Area_km2
  prediction_locations$area <- prediction_locations$Are_km2
  codes <- as.character(observations$Cod)
  dates <- as.Date(dates)

  streamflow <- demo$streamflow[
    demo$streamflow$Cod %in% codes & demo$streamflow$date %in% dates,
    c("Cod", "date", "Q")
  ]
  obs <- matrix(NA_real_, nrow = n_obs, ncol = length(dates))
  for (itime in seq_along(dates)) {
    rows <- streamflow[streamflow$date == dates[itime], ]
    rows <- rows[match(codes, rows$Cod), ]
    if (anyNA(rows$Q)) {
      stop("Missing TNDTK streamflow values in the requested fixture")
    }
    obs[, itime] <- rows$Q / observations$Area_km2
  }

  time <- as.POSIXct(dates, tz = "UTC")
  pred <- matrix(0, nrow = n_pred, ncol = length(time))

  list(
    observations = utop_stars_cube(observations, time, list(obs = obs)),
    prediction_locations = utop_stars_cube(
      prediction_locations,
      time,
      list(var1 = pred)
    ),
    streamflow = streamflow,
    codes = codes,
    time = time
  )
}
