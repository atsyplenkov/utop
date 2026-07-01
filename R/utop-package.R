#' A package providing methods for analysis and spatial interpolation of data
#' with an irregular support
#'
#' This package provides geostatistical methods for analysis and interpolation
#' of data that has an irregular support, such as as runoff characteristics or
#' population health data. The methods in this package are based on the
#' top-kriging approach suggested in Skoien et al (2006), with some
#' extensions from Gottschalk (1993).
#'
#'
#' @section Workflow:
#'
#' The work flow within the package suggests that the user is interested in a
#' prediction of a process at a series of locations where observations have not
#' been made. The example below shows a regionalization of mean annual runoff
#' in Austria.
#'
#' The easiest interface stores all variables, such as observations, prediction
#' locations, and parameters, in a [Utop] object created with
#' [utop_object()]. The `params` object below changes the default parameters so
#' that the functions use geostatistical distance instead of full
#' regularization and fit the variogram model to the variogram cloud. Most
#' user-facing functions take a `Utop` object and return an updated object.
#'
#' The data in the example below are stored as shape-files in the
#' extdata-directory of the utop-package, use the directory of your own data
#' instead. The observations consist of mean summer runoff from 138 catchments
#' in Upper Austria. The prediction locations are 863 catchments in the same
#' region. Observations and prediction locations should be stored as
#' \code{\link[sf]{sf}} polygons. Spatiotemporal data are represented as
#' vector data cubes with \code{stars}.
#'
#' \preformatted{
#' rpath = system.file("extdata", package = "utop")
#' library(sf)
#' observations = st_read(rpath, "observations")
#' prediction_locations = st_read(rpath, "predictionLocations")
#'
#' observations$obs = observations$QSUMMER_OB / observations$AREASQKM
#' params = list(g_dist_est = TRUE, g_dist_pred = TRUE, cloud = TRUE)
#' utop_obj = utop_object(observations, prediction_locations, params = params)
#' }
#'
#' The sample variogram can be added with [utop_variogram()], and a variogram
#' model can be fitted with [utop_fit_variogram()]. The fitting function calls
#' [utop_variogram()] when the object does not already contain a sample
#' variogram.
#'
#' \preformatted{
#' utop_obj = utop_variogram(utop_obj)
#' utop_obj = utop_fit_variogram(utop_obj, maxn = 2000)
#' }
#'
#' The interpolation function [utop_krige()] solves the kriging system based on
#' the computed regularized semivariances. Cross-validation can be called with
#' the argument `cv = TRUE`, either in `params` or in the call to
#' [utop_krige()].
#'
#' \preformatted{
#' utop_obj = utop_krige(utop_obj)
#' # plot predictions with sf/ggplot2
#' ggplot(utop_obj@@predictions) + aes(fill = var1.pred) + geom_sf()
#'
#' utop_obj = utop_krige(utop_obj, cv = TRUE)
#' ggplot(utop_obj@@predictions) + aes(fill = var1.var) + geom_sf()
#' }
#' @references L. Gottschalk. Interpolation of runoff applying objective
#' methods. Stochastic Hydrology and Hydraulics, 7:269-281, 1993.
#'
#' Skoien J. O., R. Merz, and G. Bloschl. Top-kriging - geostatistics on stream
#' networks. Hydrology and Earth System Sciences, 10:277-287, 2006.
#'
#' Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione,
#' A., 2014. Rtop: An R package for interpolation of data with a variable
#' spatial support, with an example from river networks. Computers &
#' Geosciences, 67.
#' @keywords spatial
#'
#' @useDynLib utop, .registration = TRUE
#' @import graphics
#' @import gstat
#' @import grDevices
#' @import stats
#' @import methods
#' @import utils
#' @importFrom units set_units
#' @import sf
#' @importFrom stars st_as_stars st_dimensions st_get_dimension_values
"_PACKAGE"
