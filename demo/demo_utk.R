# demo_utk.R
# Universal top-kriging MAF application using the local utop development
# package and altitude as a covariate.
#
# This follows demo_rtop.R, but loads the package source with
# devtools::load_all("..") and invokes universal top-kriging with a non-trivial
# right-hand side in formulaString: obs ~ Altitude. Altitude is available on the
# gauged catchment layer. The ungauged catchment layer has no altitude field, so
# this standalone demo assigns each target the nearest gauged-station altitude.
# For production use, replace this placeholder with outlet altitude,
# station altitude, or catchment-mean altitude from a DEM or trusted metadata.
#
# This file reads demo.gpkg created by demo_prep.R. Predicted/CV MAF values are
# kept in memory and are not written back to the GeoPackage.

library(sf)
devtools::load_all()
library(yardstick)
library(tidyhydro)
library(dplyr)

nselog_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  nse_vec(log(truth), log(estimate), na_rm = na_rm, ...)
}

nselog <- function(data, ...) UseMethod("nselog")
nselog <- new_numeric_metric(nselog, direction = "maximize")

nselog.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  numeric_metric_summarizer(
    name = "nselog",
    fn = nselog_vec,
    data = data,
    truth = !!rlang::enquo(truth),
    estimate = !!rlang::enquo(estimate),
    na_rm = na_rm
  )
}

input_gpkg <- "demo/demo.gpkg"

GaugedCatchments <- st_read(
  input_gpkg,
  layer = "gauged_catchments",
  quiet = TRUE
)
UngaugedCatchments <- st_read(
  input_gpkg,
  layer = "ungauged_catchments",
  quiet = TRUE
)
GaugedStations <- st_read(input_gpkg, layer = "gauged_stations", quiet = TRUE)

# Re-compute catchment areas with sf, replacing rgeos::gArea(..., byid = TRUE).
# The tutorial uses Area_km2 for gauged catchments and Are_km2 for targets.
GaugedCatchments$Area_km2 <- as.numeric(st_area(GaugedCatchments)) / 10^6
UngaugedCatchments$Are_km2 <- as.numeric(st_area(UngaugedCatchments)) / 10^6

# The ungauged catchment layer has no altitude field. Use the target outlet
# coordinates stored in the ungauged catchment attributes to assign a placeholder
# altitude from the nearest gauged station, matching demo_uk.R.
UngaugedPoints <- st_as_sf(
  st_drop_geometry(UngaugedCatchments),
  coords = c("X_3035", "Y_3035"),
  crs = st_crs(UngaugedCatchments)
)
UngaugedCatchments$Altitude <- GaugedStations$Altitud[st_nearest_feature(
  UngaugedPoints,
  GaugedStations
)]

# 6(a). Scaling relationship between MAF and drainage area for gauged catchments.
loglin.mod <- lm(log(GaugedCatchments$MAF) ~ log(GaugedCatchments$Area_km2))
c1 <- exp(loglin.mod$coefficients[1])
c2 <- loglin.mod$coefficients[2]

# Universal top-kriging prediction of MAF:
# krige obs = MAF / Area^c2 with altitude as an external drift, then rescale
# predictions by target Area^c2.
set.seed(1)
vic <- 6
GaugedCatchments$obs <- GaugedCatchments$MAF / (GaugedCatchments$Area_km2^c2)

rtop_params <- list(
  gDist = TRUE,
  rresol = 500,
  nmax = vic,
  wlim = 1,
  debug.level = 0,
  partialOverlap = TRUE
)

rtopObj.MAF <- createRtopObject(
  observations = GaugedCatchments,
  predictionLocations = UngaugedCatchments,
  formulaString = obs ~ Altitude,
  params = rtop_params
)

set.seed(21)
rtopObj.MAF <- rtopVariogram(rtopObj.MAF)
rtopObj.MAF <- rtopFitVariogram(rtopObj.MAF)
rtopObj.MAF <- checkVario(rtopObj.MAF)
rtopObj.MAF <- rtopKrige(rtopObj.MAF)

ungauged_maf <- st_drop_geometry(UngaugedCatchments) |>
  transmute(
    Locatin,
    Altitude,
    MAF_pred = rtopObj.MAF$predictions$var1.pred * (Are_km2^c2)
  )

print(ungauged_maf)

# Leave-one-out CV on MAF for gauged stations. rtopKrige(..., cv = TRUE)
# returns leave-one-out predictions for obs; convert back to dimensional MAF
# using the same area scaling as above.
rtopObj.MAF.cv <- rtopKrige(rtopObj.MAF, cv = TRUE)
cv_pred_obs <- rtopObj.MAF.cv$predictions$var1.pred

gauged_cv <- st_drop_geometry(GaugedCatchments) |>
  transmute(
    Cod = as.character(Cod),
    MAF_obs = MAF,
    MAF_pred = cv_pred_obs * (Area_km2^c2),
    MAF_resid = MAF_pred - MAF_obs
  ) |>
  filter(is.finite(MAF_obs), is.finite(MAF_pred), MAF_obs > 0, MAF_pred > 0)

maf_metrics <- metric_set(kge2012, pbias, rmse, nse, nselog)
cv_metrics <- maf_metrics(gauged_cv, truth = MAF_obs, estimate = MAF_pred)

print(cv_metrics)
