utop_spacetime_fixtures <- function(n_obs = 8, n_pred = 4, n_time = 3) {
  if (!requireNamespace("spacetime", quietly = TRUE)) {
    stop("spacetime not available")
  }

  set.seed(42)
  x <- c(0, 2, 2, 0)
  y <- c(0, 0, 2, 2)

  make_polys <- function(n, prefix) {
    polys <- lapply(seq_len(n), function(i) {
      px <- x + runif(1, 0, 20)
      py <- y + runif(1, 0, 20)
      sp::Polygons(list(sp::Polygon(cbind(px, py))), ID = paste0(prefix, i))
    })
    sp::SpatialPolygonsDataFrame(
      sp::SpatialPolygons(polys),
      data = data.frame(
        area = sapply(
          methods::slot(sp::SpatialPolygons(polys), "polygons"),
          function(i) methods::slot(i, "area")
        ),
        row.names = paste0(prefix, seq_len(n))
      )
    )
  }

  sp_obs <- make_polys(n_obs, "P")
  sp_pred <- make_polys(n_pred, "Pr")

  time <- as.POSIXct("2020-01-01") + seq_len(n_time) * 86400

  base_obs <- runif(n_obs, 1, 10)
  obs_data <- data.frame(
    obs = as.vector(sapply(base_obs, function(m) m + rnorm(n_time, 0, 0.5)))
  )
  st_obs <- spacetime::STFDF(sp_obs, time, obs_data)
  st_obs_s <- as(st_obs, "STSDF")

  # Use non-NA placeholders: as(STFDF, "STSDF") drops NA rows
  pred_data <- data.frame(var1 = rep(0, n_pred * n_time))
  st_pred <- spacetime::STFDF(sp_pred, time, pred_data)
  st_pred_s <- as(st_pred, "STSDF")

  list(observations = st_obs_s, prediction_locations = st_pred_s, time = time)
}

utop_tndtk_spacetime_fixture <- function(
  n_obs = 8,
  n_pred = 3,
  dates = as.Date("2006-01-01") + 0:4
) {
  if (!requireNamespace("spacetime", quietly = TRUE)) {
    stop("spacetime not available")
  }

  demo <- utop_demo_data()
  observations <- demo$gauged_catchments[seq_len(n_obs), ]
  prediction_locations <- demo$ungauged_catchments[seq_len(n_pred), ]
  codes <- as.character(observations$Cod)
  dates <- as.Date(dates)

  grid <- expand.grid(
    Cod = codes,
    date = dates,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  streamflow <- demo$streamflow[
    demo$streamflow$Cod %in% codes & demo$streamflow$date %in% dates,
    c("Cod", "date", "Q")
  ]
  grid <- merge(
    grid,
    streamflow,
    by = c("Cod", "date"),
    all.x = TRUE,
    sort = FALSE
  )
  grid <- grid[order(grid$date, match(grid$Cod, codes)), ]
  if (anyNA(grid$Q)) {
    stop("Missing TNDTK streamflow values in the requested fixture")
  }

  sp_obs <- as(observations, "Spatial")
  sp_obs$area <- sp_obs$Area_km2
  sp_pred <- as(prediction_locations, "Spatial")
  sp_pred$area <- sp_pred$Are_km2

  time <- as.POSIXct(dates, tz = "UTC")
  grid$obs <- grid$Q / rep(sp_obs$Area_km2, times = length(dates))
  st_obs <- spacetime::STFDF(
    sp_obs,
    time,
    data.frame(obs = grid$obs)
  )
  st_pred <- spacetime::STFDF(
    sp_pred,
    time,
    data.frame(var1 = rep(0, dim(sp_pred)[1] * length(time)))
  )

  list(
    observations = as(st_obs, "STSDF"),
    prediction_locations = as(st_pred, "STSDF"),
    streamflow = grid,
    codes = codes,
    time = time
  )
}
