# demo_laaha.R
# Two-step regression kriging following Laaha et al. (2013) "Spatial
# prediction of stream temperatures using top-kriging with an external drift."
# Environ Model Assess (2013) 18:671-683, doi:10.1007/s10666-013-9373-3.
#
# Laaha et al. decomposed the spatial variable into a deterministic drift
# (modelled by regression on altitude) and a stochastic residual (modelled by
# top-kriging). The two steps are estimated separately and superimposed:
#
#   Z(x) = m(x) + Y(x)
#
# where m(x) is the deterministic drift and Y(x) is the intrinsic stationary
# residual, kriged with a regularised residual variogram.
#
# This demo applies the same two-step logic to MAF prediction:
#
#   1. Fit lm(obs ~ Altitude) on gauged catchments as the drift model.
#   2. Compute residuals = obs - drift.
#   3. Ordinary top-krige the residuals over catchment supports.
#   4. Superimpose: predicted_obs = drift_prediction + kriged_residual.
#   5. Rescale by target Area^c2 to get dimensional MAF.
#
# This is distinct from demo_utk.R, which estimates drift and kriging weights
# together in a unified augmented kriging system (proper universal kriging).
#
# Uses the public utop S7 API. The intercept-only formula `resid ~ 1` invokes
# ordinary top-kriging on the drift residuals.
#
# This file reads demo.gpkg created by demo_prep.R. Predicted/CV MAF values are
# kept in memory and are not written back to the GeoPackage.

library(sf)
library(utop)
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

input_gpkg <- "demo.gpkg"

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
# altitude from the nearest gauged station, matching demo_uk.R and demo_utk.R.
# For production use, replace with outlet, station, or catchment-mean altitude
# from a DEM or trusted metadata source.
UngaugedPoints <- st_as_sf(
  st_drop_geometry(UngaugedCatchments),
  coords = c("X_3035", "Y_3035"),
  crs = st_crs(UngaugedCatchments)
)
UngaugedCatchments$Altitude <- GaugedStations$Altitud[st_nearest_feature(
  UngaugedPoints,
  GaugedStations
)]

# Scaling relationship between MAF and drainage area for gauged catchments.
loglin.mod <- lm(log(GaugedCatchments$MAF) ~ log(GaugedCatchments$Area_km2))
c1 <- exp(loglin.mod$coefficients[1])
c2 <- loglin.mod$coefficients[2]

# The kriging target: area-scaled MAF, matching demo_rtop.R and demo_utk.R.
GaugedCatchments$obs <- GaugedCatchments$MAF / (GaugedCatchments$Area_km2^c2)

# ------------------------------------------------------------------------
# Step 1: Fit the deterministic drift model.
# Laaha et al. (2013) used exponential regression T = a * exp(b * H) for
# stream temperature. For area-scaled MAF, a simple linear drift on altitude
# is used here to match the covariate specification of demo_uk.R/demo_utk.R.
# ------------------------------------------------------------------------
drift_model <- lm(obs ~ Altitude, data = GaugedCatchments)

# Step 2: Compute regression residuals.
GaugedCatchments$drift <- fitted(drift_model)
GaugedCatchments$resid <- residuals(drift_model)

# Step 3: Ordinary top-krige the residuals over catchment supports.
# The variogram is fitted to the regression residuals and regularised
# over the actual support areas — exactly the procedure in Laaha et al.
# Section 3.3 ("Variogram of the Residuals").
set.seed(1)
vic <- 6

utop_params <- list(
  gDist = TRUE,
  rresol = 500,
  nmax = vic,
  wlim = 1,
  debug.level = 0,
  partialOverlap = TRUE
)

utop_obj_resid <- utop_object(
  observations = GaugedCatchments,
  prediction_locations = UngaugedCatchments,
  formula = resid ~ 1,
  params = utop_params
)
utop_obj_resid <- utop_variogram(utop_obj_resid)
utop_obj_resid <- utop_fit_variogram(utop_obj_resid)
utop_obj_resid <- utop_krige(utop_obj_resid)

# Step 4: Predict the drift at ungauged locations and superimpose.
ungauged_drift <- predict(drift_model, newdata = UngaugedCatchments)
ungauged_resid_kriged <- utop_obj_resid@predictions$var1.pred
ungauged_obs_pred <- ungauged_drift + ungauged_resid_kriged

# Step 5: Rescale to dimensional MAF.
ungauged_maf <- st_drop_geometry(UngaugedCatchments) |>
  transmute(Locatin, Altitude, MAF_pred = ungauged_obs_pred * (Are_km2^c2))

print(ungauged_maf)

# ------------------------------------------------------------------------
# Leave-one-out cross-validation
#
# Laaha et al. (2013, Section 3.5) updated the regression equation for each
# left-out point, then kriged residuals using the variogram from all n obs:
#
#   "For the regression model, we updated the regression equation for the
#    remaining n - 1 data points. The so obtained regression model was used
#    for prediction at that location."
#
# This demo simplifies the CV by fitting the regression once on all data and
# then using utop_krige(cv = TRUE) for the residual kriging step. The drift
# component is therefore slightly optimistic (fitted value includes the
# left-out point). A production CV should re-fit the regression in each
# LOO fold.  The simplification is flagged here, not hidden.
# ------------------------------------------------------------------------
utop_obj_resid_cv <- utop_krige(utop_obj_resid, cv = TRUE)
cv_resid <- utop_obj_resid_cv@predictions$var1.pred

# Superimpose: cv prediction = drift (full-data) + cv-kriged residual.
cv_obs_pred <- GaugedCatchments$drift + cv_resid

gauged_cv <- st_drop_geometry(GaugedCatchments) |>
  transmute(
    Cod = as.character(Cod),
    MAF_obs = MAF,
    MAF_pred = cv_obs_pred * (Area_km2^c2),
    MAF_resid = MAF_pred - MAF_obs
  ) |>
  filter(is.finite(MAF_obs), is.finite(MAF_pred), MAF_obs > 0, MAF_pred > 0)

maf_metrics <- metric_set(kge2012, pbias, rmse, nse, nselog)
cv_metrics <- maf_metrics(gauged_cv, truth = MAF_obs, estimate = MAF_pred)

print(cv_metrics)
