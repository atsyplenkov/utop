# Demo material

The current demo workflow is based on the TNDTK example by Simone Persiano and collaborators:

- GitHub repository: <https://github.com/SimonePersiano/TNDTK>
- Original tutorial <https://rpubs.com/simone_persiano/927821>

The demo uses the Tyrol / South Tyrol subset distributed with the TNDTK tutorial: 27 gauged catchments, 27 gauging stations, 3 ungauged target catchments, and daily streamflow. Daily streamflow were used to estimate MAF. The `demo.gpkg` is a complete copy of the shapefiles used by Persiano et al. with the only difference of hardcoded MAF values. 

### `demo_rtop.R`

Runs MAF top-kriging using [`rtop`](https://cran.r-project.org/package=rtop). This follows the MAF branch of TNDTK Tutorial:

1. Fit the log-linear scaling relationship between observed MAF and catchment area.
2. Krige `obs = MAF / Area^c2` with catchment polygons as supports.
3. Rescale predictions by target catchment area.
4. Run leave-one-out CV for gauged catchments.
5. Print metrics with `yardstick::metric_set()` and `tidyhydro` metrics.

### `demo_utk.R`

Runs universal top-kriging with the public `utop` S7 API and altitude as an external drift term:

```r
obs ~ Altitude
```

This is the same area-scaled MAF workflow as `demo_rtop.R`: `obs = MAF / Area^c2` is kriged over catchment supports and predictions are rescaled by target catchment area. The `ungauged_catchments` layer has no altitude field, so the script assigns each target the nearest gauged-station altitude as a placeholder.

### `demo_ok.R`

Runs an ordinary-kriging baseline using [`automap`](https://cran.r-project.org/package=automap) on gauging station **points**.

The hydrological variable follows Farmer (2016):

```r
z = log(MAF / drainage_area)
```

Predictions are back-transformed and multiplied by target drainage area.

### `demo_uk.R`

Runs universal kriging using [`automap`](https://cran.r-project.org/package=automap) on gauging station **points**, with altitude as a covariate/drift term:

```r
log_unit_maf ~ Altitud
```

The `ungauged_catchments` layer has no altitude field. For this standalone demo, the script assigns each target the nearest gauged-station altitude. This is a placeholder. For production use, replace it with outlet altitude, station altitude, or catchment-mean altitude from a DEM or trusted metadata source.

### `demo_laaha.R`

Runs the two-step regression kriging approach of Laaha et al. (2013) — **not** the unified universal top-kriging of `demo_utk.R`. The procedure follows Laaha et al. Section 2.2 ("Extending Top-Kriging with an External Drift Function"):

1. Fit `lm(obs ~ Altitude)` on gauged catchments as the deterministic drift model.
2. Compute residuals `resid = obs - drift`.
3. Ordinary top-krige the residuals over catchment supports.
4. Superimpose: `predicted_obs = drift + kriged_residual`.
5. Rescale by target `Area^c2` to get dimensional MAF.

Uses the public `utop` S7 API for the top-kriging step with `formula = resid ~ 1`. The drift and kriging weights are estimated **separately** — the regression first, then the residual variogram and kriging. In `demo_utk.R`, by contrast, the trend and kriging weights are solved **together** in a single augmented system.

**CV simplification note.** Laaha et al. (2013, Section 3.5) re-fit the regression for each left-out point. This demo fits the regression once on all data and then uses `utop_krige(cv = TRUE)` on the residuals. The drift component in CV is therefore slightly optimistic. A production comparison should re-fit the regression in each LOO fold; the simplification is documented, not hidden.

## Leave-one-out CV comparison on gauged MAF

| Method | KGE2012 | PBIAS | RMSE | NSE | NSElog |
|---|---:|---:|---:|---:|---:|
| Top-Kriging | 0.963 | 2.45 | 0.697 | 0.978 | 0.945 |
| Universal Top-Kriging | 0.971 | 2.04 | 0.673 | 0.979 | 0.940 |
| Laaha Two-Step TK+ED | 0.957 | 0.333 | 0.833 | 0.968 | 0.937 |
| Ordinary Kriging | 0.960 | 0.162 | 0.916 | 0.961 | 0.929 |
| Universal Kriging | 0.978 | -1.60 | 0.816 | 0.969 | 0.955 |

## Ungauged MAF predictions (approximate because ungauged target altitude is unknown)

| Target | Top-Kriging | Universal Top-Kriging | Laaha Two-Step TK+ED | Ordinary Kriging | Universal Kriging |
|---|---:|---:|---:|---:|---:|
| Ahr_3 | 10.436 | 10.442 | 9.977 | 11.454 | 10.766 |
| Gader_1 | 6.178 | 6.303 | 5.481 | 5.789 | 5.512 |
| Isel_4 | 2.226 | 2.303 | 2.419 | 2.153 | 2.245 |


## References

- Castellarin, A., Persiano, S., Pugliese, A., Aloe, A., Skøien, J. O., and Pistocchi, A. (2018). Prediction of streamflow regimes over large geographical areas: interpolated flow-duration curves for the Danube Region. *Hydrological Sciences Journal*, 63(6), 845-861. <https://doi.org/10.1080/02626667.2018.1445855>
- Persiano, S., Pugliese, A., Aloe, A., Skøien, J. O., Castellarin, A., and Pistocchi, A. (2022). Streamflow data availability in Europe: a detailed dataset of interpolated flow-duration curves. *Earth System Science Data Discussions*. <https://doi.org/10.5194/essd-2021-469>
- Pugliese, A., Castellarin, A., and Brath, A. (2014). Geostatistical prediction of flow-duration curves in an index-flow framework. *Hydrology and Earth System Sciences*, 18, 3801-3816. <https://doi.org/10.5194/hess-18-3801-2014>
- Pugliese, A., Farmer, W. H., Castellarin, A., Archfield, S. A., and Vogel, R. M. (2016). Regional flow duration curves: geostatistical techniques versus multivariate regression. *Advances in Water Resources*, 96, 11-22. <https://doi.org/10.1016/j.advwatres.2016.06.008>
- Pugliese, A., Persiano, S., Bagli, S., Mazzoli, P., Parajka, J., Arheimer, B., Capell, R., Montanari, A., Blöschl, G., and Castellarin, A. (2018). A geostatistical data-assimilation technique for enhancing macro-scale rainfall-runoff simulations. *Hydrology and Earth System Sciences*, 22, 4633-4648. <https://doi.org/10.5194/hess-22-4633-2018>
- Skøien, J. O., Merz, R., and Blöschl, G. (2006). Top-kriging: geostatistics on stream networks. *Hydrology and Earth System Sciences*, 10, 277-287. <https://doi.org/10.5194/hess-10-277-2006>
- Farmer, W. H. (2016). Ordinary kriging as a tool to estimate historical daily streamflow records. *Hydrology and Earth System Sciences*, 20, 2721-2735. <https://doi.org/10.5194/hess-20-2721-2016>
- Laaha, G., Skøien, J. O., Nobilis, F., and Blöschl, G. (2013). Spatial prediction of stream temperatures using top-kriging with an external drift. *Environmental Modeling & Assessment*, 18, 671-683. <https://doi.org/10.1007/s10666-013-9373-3>
