# this is needed to make R aware we're introducing new S3 methods:

#' create variogram for data with spatial support
#'
#' rtopVariogram will create binned variogram or cloud variogram of data with
#' an areal support.
#'
#'
#' @param object object of class \code{rtop} (see \link{utop-package}) or a
#' \cr \code{\link[sp:SpatialPolygons]{SpatialPolygonsDataFrame}} or
#' \code{\link[sp:SpatialPoints]{SpatialPointsDataFrame}} with information
#' about observations. If \cr \code{object} is a \cr
#' \code{\link[sp:SpatialPoints]{SpatialPointsDataFrame}}, it must have a
#' column with name \code{area}.
#' @param formulaString formula that defines the dependent variable as a linear
#' model of independent variables; suppose the dependent variable has name
#' \code{z}, for ordinary and simple kriging use the formula \code{z~1}; for
#' universal kriging, suppose \code{z} is linearly dependent on \code{x} and
#' \code{y}, use the formula \code{z~x+y}. The formulaString defaults to
#' \code{"value~1"} if \code{value} is a part of the data set.  If not, the
#' first column of the data set is used.
#' @param params a set of parameters, used to modify the default parameters for
#' the \code{rtop} package, set in \code{\link{getRtopParams}}.
#' @param cloud logical; if TRUE, calculate the semivariogram cloud, can be
#' used to overrule the cloud parameter in params.
#' @param abins possibility to set areal bins (not yet implemented)
#' @param dbins possibility to set distance bins (not yet implemented)
#' @param data.table an option to use data.table internally for the variogram
#' computation for \code{\link[spacetime]{STSDF}}-objects
#' @param discPoints list of discretisation points of the observation areas
#' from \code{\link{rtopDisc}}, only used for block-averaging of
#' coordinate-based trend basis functions when the formula describes
#' universal kriging and \code{ukTrendSupport = "block"}; computed internally
#' when needed and not supplied
#' @param ... parameters to other functions called, e.g. gstat's
#' \code{\link[gstat]{variogram}}-function and to
#' \code{\link{rtopVariogram.SpatialPointsDataFrame}} when the method is called
#' with an object of a different class
#' @return The function creates a variogram, either of type
#' \code{rtopVariogram} or \code{rtopVariogramCloud}. This variogram is based
#' on the \code{\link[gstat]{variogram}} function from gstat, but has
#' additional information about the spatial size or length of the observations.
#' An rtop-object with the variogram added is returned if the function is
#' called with an rtop-object as argument.
#'
#' For spatio-temporal objects (\code{\link[spacetime]{STSDF}}), the variogram
#' is the spatially variogram, averaged for all time steps. There is a
#' possibility to use data.table internally in this function, which can improve
#' computation time for some cases.
#' @note The variogram cloud is similar to the variogram cloud from
#' \code{\link[gstat]{gstat}}, with the area/length added to the resulting
#' data.frame. The binned variogram is also based on the area or length, in
#' addition to the distance between observations. The bins equally distanced in
#' the log10-space of the distances and areas (lengths). The size of the bins
#' is decided from the parameters \code{amul} and \code{dmul}, defining the
#' number of bins per order of magnitude (1:10, 10:100, and so on).
#'
#' The distances between areas are in this function based on the centre of
#' gravity.
#' @author Jon Olav Skoien
#' @references Skoien J. O., R. Merz, and G. Bloschl. Top-kriging -
#' geostatistics on stream networks. Hydrology and Earth System Sciences,
#' 10:277-287, 2006.
#'
#' Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione,
#' A., 2014. Rtop: An R package for interpolation of data with a variable
#' spatial support, with an example from river networks. Computers &
#' Geosciences, 67.
#' @keywords spatial
#' @examples
#'
#' \dontrun{
#' library(sf)
#' rpath <- system.file("extdata",package="utop")
#' observations <- sf::st_read(rpath,"observations")
#' # Create a column with the specific runoff:
#' observations$obs <- observations$QSUMMER_OB/observations$AREASQKM
#'
#' vario <- rtopVariogram(observations, cloud = TRUE)
#' }
#'
#' @export
rtopVariogram <- function(object, ...) UseMethod("rtopVariogram")


#' Fit variogram model to sample variogram of data with spatial support
#'
#' rtopFitVariogram will fit a variogram model to the estimated binned
#' variogram or cloud variogram of data with an areal support.
#'
#'
#' @param object object of class \code{rtopVariogram} or
#' \code{rtopVariogramCloud}, or an object with class \code{rtop} that includes
#' the sample variograms.
#'
#' The object can also be of class \code{\link[sp]{SpatialPolygonsDataFrame}}
#' or \cr \code{\link[sp]{SpatialPointsDataFrame}} with observations. If
#' \code{object} is a \cr \code{\link[sp]{SpatialPointsDataFrame}}, it must
#' have a column with name \code{area}.
#' @param observations the observations, passed as a Spatial*DataFrame object,
#' if object is an \cr \code{rtopVariogram} or \code{rtopVariogramCloud}
#' @param params a set of parameters, used to modify the default parameters for
#' the \code{rtop} package, set in \code{\link{getRtopParams}}. The argument
#' params can also be used for the other methods, through the ...-argument.
#' @param dists either a matrix with geostatistical distances (created by a
#' call to the function \code{\link{gDist}} or a list with the areas
#' discretized (from a call to \code{\link{rtopDisc}}.
#' @param mr logical; defining whether the function should return a list with
#' discretized elements and geostatistical distances, even if it was not called
#' with an rtop-object as argument.
#' @param aOver a matrix with the overlapping areas of the observations, used
#' for computation of the nugget effect.  It will normally be recomputed by the
#' function if it is NULL and necessary
#' @param iprint print flag that is passed to \code{sceua::sceua}
#' @param ... Other parameters to functions called from
#' \code{rtopFitVarigoram}. For the three first methods of the function,
#' \code{...} can also include parameters to the last two methods.
#' @return The function creates an object with the fitted variogram Model
#' (\code{variogramModel}) and a \cr \code{\link{data.frame}} (\code{varFit})
#' with the differences between the sample semivariances and the regularized
#' semivariances. If \code{mr} = TRUE, the function also returns other objects
#' (discretized elements and geostatistical distances, if created) as a part of
#' the returned object. If the function is called with an rtop-object as
#' argument, it will return an rtop-object with \code{variogramModel} and
#' \code{varFit} added to the object, in addition to other objects created.
#' @note There are several options for fitting of the variogramModel, where the
#' parameters can be set in \code{params}, which is a list of parameters for
#' modification of the default parameters of the utop-package given in a call
#' to \code{\link{getRtopParams}}. The first choice is between individual
#' fitting and binned fitting. This is based on the type of variogram
#' submitted, individual fitting is done if a cloud variogram (of class
#' \code{rtopVariogramCloud}) is passed as argument, and binned fitting if the
#' submitted variogram is of class \code{rtopVariogram}. If the function is
#' called with an object of class \code{rtop}, having both \code{variogram} and
#' \code{variogramCloud} among its arguments, the variogram model is fitted to
#' the variogram which is consistent with the parameter \code{cloud}.
#' @author Jon Olav Skoien
#' @references Skoien J. O., R. Merz, and G. Bloschl. Top-kriging -
#' geostatistics on stream networks. Hydrology and Earth System Sciences,
#' 10:277-287, 2006.
#'
#' Skoien, J. O. and G. Bloschl. Spatio-Temporal Top-Kriging of Runoff Time
#' Series. Water Resources Research 43:W09419, 2007.
#'
#' Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione,
#' A., 2014. Rtop: An R package for interpolation of data with a variable
#' spatial support, with an example from river networks. Computers &
#' Geosciences, 67.
#' @keywords spatial
#' @examples
#'
#' \donttest{
#' rpath <- system.file("extdata",package="utop")
#' library(sf)
#' observations <- sf::st_read(rpath, "observations")
#' predictionLocations <- sf::st_read(rpath,"predictionLocations")
#'
#' # Create a column with the specific runoff:
#' observations$obs <- observations$QSUMMER_OB/observations$AREASQKM
#'
#' # Setting some parameters
#' params <- list(gDist = TRUE, cloud = FALSE)
#' # Create a column with the specific runoff:
#' observations$obs <- observations$QSUMMER_OB/observations$AREASQKM
#' # Build an object
#' rtopObj <- createRtopObject(observations,predictionLocations,
#'            params = params)
#' # Fit a variogram (function also creates it)
#' rtopObj <- rtopFitVariogram(rtopObj)
#' rtopObj$variogramModel
#' }
#'
#' @export
rtopFitVariogram <- function(object, ...) UseMethod("rtopFitVariogram")
#estimateParameters <- function(object, ...) UseMethod("estimateParameters")
#spatialPredict      <- function(object, ...) UseMethod("spatialPredict")

#' Plot variogram fitted to data with support
#'
#' The function will create diagnostic plots for analysis of the variograms
#' fitted to sample variograms of data with support
#'
#'
#' @param object either: object of class \code{rtop} (see
#' \code{\link{utop-package}}), or an object of type \cr \code{rtopVariogram}
#' @param acor unit correction factor in the key, e.g. to see numbers more
#' easily interpretable for large areas. As an example, ucor = 0.000001 when
#' area is given in square meters and should rather be shown as square
#' kilometers.  Note that this parameter also changes the value of the nugget
#' to the new unit.
#' @param log text variable for log-plots, default to log-log \code{"xy"}, can
#' otherwise be set to \code{"x"}, \code{"y"} or \code{""}
#' @param cloud logical; whether to look at the cloud variogram instead of the
#' binned variogram
#' @param gDist logical; whether to use ghosh-distance for semivariogram
#' regularization instead of full integration of the semivariogram
#' @param sampleVariogram a sample variogram of the data
#' @param observations a set of observations
#' @param areas either an array of areas that should be used as examples, or
#' the number of areas per order of magnitude (similar to the parameter
#' \code{amul}; see \code{\link{getRtopParams}}. amul from \code{rtopObj} or
#' from the default parameter set will be used if not defined here.
#' @param dists either an array of distances that should be used as examples,
#' or the number of distances per order of magnitude(similar to the parameter
#' \code{amul}; see \code{\link{getRtopParams}}. amul from \code{rtopObj} or
#' from the default parameter set will be used if not defined here.
#' @param acomp either a matrix with the area bins that should be visualized,
#' or a number giving the number of pairs to show. If a sample variogram is
#' given, the \code{acomp} pairs with highest number of pairs will be used
#' @param curveSmooth logical or numerical; describing whether the curves in
#' the last plot should be smoothed or not. If numeric, it gives the degrees of
#' freedom (df) for the splines used for smoothing. See also
#' \code{\link[stats]{smooth.spline}}
#' @param params list of parameters to modify the default parameters of rtopObj
#' or the default parameters found from \code{\link{getRtopParams}}
#' @param compVars a list of variograms of \code{gstat}-type for comparison,
#' see \code{\link[gstat:vgm]{vgm}}. The names of the variograms in the list
#' will be used in the key.
#' @param legx x-coordinate of the legend for fine-tuning of position, see
#' x-argument of \cr \code{\link[graphics]{legend}}
#' @param legy y-coordinate of the legend for fine-tuning of position, see
#' y-argument of \cr \code{\link[graphics]{legend}}
#' @param plotNugg logical; whether the nugget effect should be added to the
#' plot or not
#' @param ... arguments to lower level functions
#' @return The function gives diagnostic plots for the fitted variograms, where
#' the regularized variograms are shown together with the sample variograms and
#' possibly also user defined variograms. In addition, if an rtopObject is
#' submitted, the function will also give plots of the relationship between
#' variance and area size and a scatter plot of the fit of the observed and
#' regularized variogram values. The sizes of the dots are relative to the
#' number of pairs in each group.
#' @author Jon Olav Skoien
#' @references Skoien J. O., R. Merz, and G. Bloschl. Top-kriging -
#' geostatistics on stream networks. Hydrology and Earth System Sciences,
#' 10:277-287, 2006.
#'
#' Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione,
#' A., 2014. Rtop: An R package for interpolation of data with a variable
#' spatial support, with an example from river networks. Computers &
#' Geosciences, 67.
#' @keywords spatial
#' @examples
#'
#' \donttest{
#' library(gstat)
#' rpath <- system.file("extdata",package="utop")
#' library(sf)
#' observations <- sf::st_read(rpath, "observations")
#' predictionLocations <- sf::st_read(rpath,"predictionLocations")
#'
#' # Create a column with the specific runoff:
#' observations$obs <- observations$QSUMMER_OB/observations$AREASQKM
#' params <- list(cloud = TRUE, gDist = TRUE)
#' rtopObj <- createRtopObject(observations, predictionLocations,
#'                            params = params)
#'
#' # Fit a variogram (function also creates it)
#' rtopObj <- rtopFitVariogram(rtopObj, maxn = 2000)
#' checkVario(rtopObj,
#'     compVars = list(first = gstat::vgm(5e-6, "Sph", 30000,5e-8),
#'                    second = gstat::vgm(2e-6, "Sph", 30000,5e-8)))
#'
#' rtopObj <- checkVario(rtopObj, acor = 0.000001,
#'           acomp = data.frame(acl1 = c(2,2,2,2,3,3,3,4,4),
#'           acl2 = c(2,3,4,5,3,4,5,4,5)))
#' rtopObj <- checkVario(rtopObj, cloud = TRUE, identify = TRUE,
#'           acor = 0.000001)
#' }
#'
#' @export
checkVario <- function(object, ...) UseMethod("checkVario")

#' calculate geostatistical distances between areas
#'
#' Calculate geostatistical distances (Ghosh-distances) between areas
#'
#'
#' @param object object of class \code{\link[sp]{SpatialPolygons}} or
#' \code{\link{rtopDisc}}; or object of class \code{rtop} with such boundaries
#' and/or discretized elements (the individual areas)
#' @param params a set of parameters, used to modify the default parameters for
#' the \code{rtop} package, set in \code{\link{getRtopParams}}. The argument
#' params can also be used for the other methods, through the ...-argument.
#' @param object2 an object of same type as \code{object}, except for
#' \code{rtop}; for calculation of geostatistical distances also between the
#' elements in the two different objects
#' @param diag logical; if TRUE only calculate the geostatistical distances
#' between each element and itself, only when the objects are lists of
#' discretized areas and object2 = object or object2 = NULL
#' @param debug.level debug.level = 0 will suppress output from the call to
#' varMat, done for calculation of the geostatistical distances
#' @param ... other parameters, for \code{gDist.list} when calling one of the
#' other methods, or for \code{\link{varMat}}, in which the calculations take
#' place
#' @return If called with one list of discretized elements, a matrix with the
#' geostatistical distances between the elements within the list. If called
#' with two lists of discretized elements, a matrix with the geostatistical
#' distances between the elements in the two lists. If called with \code{diag =
#' TRUE}, the function returns an array of the geostatistical distance within
#' each of the elements in the list.
#'
#' If called with one \code{\link[sp]{SpatialPolygons}} or
#' \code{\link[sp]{SpatialPolygonsDataFrame}} or the function returns a list
#' with one matrix with geostatistical distances between the elements of the
#' object. If called with two objects, the list will also containt a matrix of
#' the geostatistical distances between the elements of the two objects, and an
#' array of the geostatistical distances within the elements of the second
#' object.
#'
#' If called with an rtop-object, the function will return the object, amended
#' with the list above.
#' @note The geostatistical distance can be seen as the average distance
#' between points in two elements, or the average distance within points in a
#' single element. The distance measure is also sometimes referred to as
#' Ghosh-distance, from Ghosh (1951) who found analytical expressions for these
#' distances between blocks with regular geometry.
#'
#' The use of geostatistical distances within \code{rtop} is based on an idea
#' from Gottschalk (1993), who suggested to replace the traditional
#' regularization of variograms within block-kriging (as done in the original
#' top-kriging application of Skoien et al (2006)) with covariances of the
#' geostatistical distance. The covariance between two areas can then be found
#' as \code{C(a1,a2) = cov(gd)} where \code{gd} is the geostatistical distance
#' between the two areas \code{a1} and \code{a2}, instead of an integration of
#' the covariance function between the two areas.
#'
#' \code{rtop} is based on semivariograms instead of covariances, and the
#' semivariogram value between the two areas can be found as \code{gamma(a1,a2)
#' = g(gd) - 0.5 (g(gd1) + g(gd2))} where \code{g} is a semivariogram valid for
#' point support, \code{gd1)} and \code{gd2} are the geostatistical distances
#' within each of the two areas.
#' @author Jon Olav Skoien
#' @references
#'
#' Ghosh, B. 1951. Random distances within a rectangle and between two
#' rectangles. Bull. Calcutta Math. Soc., 43, 17-24.
#'
#' Gottschalk, L. 1993. Correlation and covariance of runoff. Stochastic
#' Hydrology and Hydraulics, 7, 85-101.
#'
#' Skoien, J. O., R. Merz, and G. Bloschl. 2006. Top-kriging - geostatistics on
#' stream networks. Hydrology and Earth System Sciences, 10, 277-287.
#'
#' Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione,
#' A., 2014. Rtop: An R package for interpolation of data with a variable
#' spatial support, with an example from river networks. Computers &
#' Geosciences, 67.
#' @keywords spatial
#' @examples
#'
#' \donttest{
#' rpath <- system.file("extdata",package="utop")
#' library(sf)
#' observations <- sf::st_read(rpath, "observations")
#' gDist <- gDist(observations)
#' }
#'
#' @export
gDist <- function(object, ...) UseMethod("gDist")


#' Discretize areas
#'
#' \code{rtopDisc} will discretize an area for regularization or calculation of
#' Ghosh-distance
#'
#' There are different options for discretizing the objects. When the areas
#' from the bins are discretized, the options are \code{random} or
#' \code{regular} sampling, \code{regular} sampling is the default.
#'
#' For the real areas, regular sampling appears to have computational
#' advantages compared with random sampling. In addition to the traditional
#' regular sampling, \code{rtop} also offers a third type of sampling which
#' assures that the same discretization points are used for overlapping areas.
#'
#' Starting with a coarse grid covering the region of interest, this will for a
#' certain support be refined till a requested minimum number of points from
#' the grid is within the support.  In this way, for areal supports, the number
#' of points in the area with the largest number of points will be
#' approximately four times the requested minimum number of points. This
#' methods also assure that points used to discretize a large support will be
#' reused when discretizing smaller supports within the large one, e.g.
#' subcatchments within larger catchments.
#'
#' @param object object of class \code{\link[sp]{SpatialPolygons}} or
#' \code{\link[sp:SpatialPolygons]{SpatialPolygonsDataFrame}} or
#' \code{rtopVariogram}, or an object with class \code{rtop} that includes one
#' of the above
#' @param bb boundary box, usually modified to be the common boundary box for
#' two spatial object
#' @param params possibility to pass parameters to modify the default
#' parameters for the \code{rtop} package, set in \code{\link{getRtopParams}}.
#' Typical parameters to modify for this function are:
#' - rresol = 100; minimum number of discretization points in areas
#' - hresol = 5; number of discretization points in one direction for areas
#'   in binned variograms
#' - hstype = "regular"; sampling type for binned variograms
#' - rstype = "rtop"; sampling type for real areas
#' @param ... Possibility to pass individual parameters
#' @return The function returns a list of discretized areas, or if called with
#' an rtop-object as argument, the object with lists of discretizations of the
#' observations and prediction locations (if part of the object). If the
#' function is called with an rtopVariogram (usually this is an internal call),
#' the list contains discretized pairs of hypothetical objects from each bin of
#' the semivariogram with a centre-to-centre distance equal to the average
#' distance between the objects in a certain bin.
#' @author Jon Olav Skoien
#' @seealso \code{\link{rtopVariogram}}
#' @references Skoien J. O., R. Merz, and G. Bloschl. Top-kriging -
#' geostatistics on stream networks. Hydrology and Earth System Sciences,
#' 10:277-287, 2006.
#'
#' Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione,
#' A., 2014. Rtop: An R package for interpolation of data with a variable
#' spatial support, with an example from river networks. Computers &
#' Geosciences, 67.
#' @keywords spatial
#' @export
rtopDisc <- function(object, ...) UseMethod("rtopDisc")

#' create a semivariogram matrix between a set of locations, or semivariogram
#' matrices between and within two sets of locations
#'
#' varMat will create a semivariogram matrix between all the supports in a set
#' of locations (observations or prediction locations) or semivariogram
#' matrices between all the supports in one or two sets of locations, and also
#' between them.
#'
#'
#' @param object either: 1) an object of class \code{rtop} (see
#' \link{utop-package}) or 2) a \cr
#' \code{\link{matrix}} with geostatistical distances (see \code{\link{gDist}}
#' or 4) a \code{\link{list}} with discretized supports
#' @param varMatUpdate logical; if TRUE, also existing variance matrices will
#' be recomputed, if FALSE, only missing variance matrices will be computed
#' @param fullPred logical; whether to create the full covariance matrix also
#' for the predictions, mainly used for simulations
#' @param params a set of parameters, used to modify the default parameters for
#' the \code{rtop} package, set in \code{\link{getRtopParams}}.
#' @param object2 if \code{object} is not an object of class \code{rtop}; an
#' object of the same class as \code{object} with a possible second set of
#' locations with support
#' @param variogramModel variogramModel to be used in calculation of the
#' semivariogram matrix (matrices)
#' @param ... typical parameters to modify from the default parameters of the
#' utop-package (or modifications of the previously set parameters for the
#' \code{rtop}-object), see also \code{\link{getRtopParams}}. These can also be
#' passed in a list named params, as for the rtop-method. Typical parameters to
#' modify for this function: \describe{ \item{rresol = 100}{miminum number of
#' discretization points, in call to \code{\link{rtopDisc}} if necessary}
#' \item{rstype = "rtop"}{sampling type from areas, in call to
#' \code{\link{rtopDisc}} if necessary} \item{gDistPred = FALSE}{use
#' geostatistical distance for semivariogram matrices} \item{gDist}{parameter
#' to set jointly \code{gDistEst = gDistPred = gDist}} }
#' @param overlapObs matrix with observations that overlap each other
#' @param overlapPredObs matrix with \code{observations} and
#' \code{predictionLocations} that overlap each other
#' @param coor1 coordinates of centroids of \code{object}
#' @param coor2 coordinates of centre-of-gravity of \code{object2}
#' @param maxdist maximum distance between areas for inclusion in semivariogrma
#' matrix
#' @param diag logical; if TRUE only the semivariogram values along the
#' diagonal will be calculated, typical for semivariogram matrix of prediction
#' locations
#' @param sub1 semivariogram array for subtraction of inner variances of areas
#' @param sub2 semivariogram array for subtraction of inner variances of areas
#' @param debug.level debug.level >= 1 will give output for every element
#' @return The lower level versions of the function calculates a semivariogram
#' matrix between locations in \code{object} or between the locations in
#' \code{object} and the locations in \code{object2}. The method for object of
#' type \code{rtop} calculates semivariogram matrices between observation
#' locations, between prediction locations, and between observation locations
#' and prediction locations, and adds these to \code{object}.
#' @note The argument \code{varMatUpdate} is typically used to avoid repeated
#' computations of the same variance matrices. Default is FALSE, which will
#' avoid recomputation of the variance matrix for the observations if the
#' procedure is cross-validation before interpolation. Should be set to TRUE if
#' the variogram Model has been changed, or if observation and/or prediction
#' locations have been changed.
#'
#' If an \code{rtop}-object contains observations and/or predictionLocations of
#' type \code{\link[spacetime]{STSDF}}, the covariance matrix is computed based
#' on the spatial properties of the object.
#' @author Jon Olav Skoien
#' @seealso \code{\link{gDist}}
#' @references Skoien J. O., R. Merz, and G. Bloschl. Top-kriging -
#' geostatistics on stream networks. Hydrology and Earth System Sciences,
#' 10:277-287, 2006.
#'
#' Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione,
#' A., 2014. Rtop: An R package for interpolation of data with a variable
#' spatial support, with an example from river networks. Computers &
#' Geosciences, 67.
#' @keywords spatial
#' @examples
#'
#' \dontrun{
#' library(sf)
#' rpath <- system.file("extdata",package="utop")
#' observations <- sf::st_read(rpath,"observations")
#' vmod <- list(model = "Ex1", params = c(0.00001,0.007,350000,0.9,1000))
#' vm <- varMat(observations, variogramModel = vmod)
#' }
#'
#' @export
varMat <- function(object, ...) UseMethod("varMat")


#' Spatial interpolation of data with spatial support
#'
#' rtopKrige perform spatial interpolation or cross validation of data with
#' areal support.
#'
#' This function is the interpolation routine of the utop-package. The
#' simplest way of calling the function is with an rtop-object that contains
#' the fitted variogram model and all the other necessary data (see
#' \code{\link{createRtopObject}} or \code{\link{utop-package}}).
#'
#' The function will, if called with covariance matrices between observations
#' and between observations and prediction locations, use these for the
#' interpolation. If the function is called without these matrices,
#' \code{\link{varMat}} will be called to create them. These matrices can
#' therefore be reused if necessary, an advantage as it is computationally
#' expensive to create them.
#'
#' The interpolation that takes part within \code{rtopKrige.default} is based
#' on the semivariance matrices between observations and between observations
#' and prediction locations. It is therefore possible to use this function also
#' to interpolate data where the matrices have been created in other ways, e.g.
#' based on distances in physiographical space or distances along a stream.
#'
#' The function returns the weights rather than the predictions if \code{wret =
#' TRUE}. This is useful for batch processing of time series, e.g. once the
#' weights are created, they can be used to compute the interpolated values for
#' each time step.
#'
#' rtop is able to take some advantage of multiple CPUs, which can be invoked
#' with the parameter \code{nclus}. When it gets a number larger than one,
#' \code{rtopKrige} will start a cluster with \code{nclus} workers, if the
#' \code{\link{parallel}}-package has been installed.
#'
#' The parameter \code{singularSolve} can be used when some areas are almost
#' completely overlapping. In this case, the discretization of them might be
#' equal, and the covariances to other areas will also be equal. The kriging
#' matrix will in this case be singular. When \code{singularSolve = TRUE},
#' \code{rtopKrige} will remove one of the neighbours, and instead work with
#' the mean of the two observations. An overview of removed neighbours can be
#' seen in the resulting object, under the name \code{removed}.
#'
#' Kriging of time series is possible when \code{observations} and
#' \code{predictionLocations} are spatiotemporal objects of type
#' \code{\link[spacetime]{STSDF}}. The interpolation is still spatial, in the
#' sense that the regularisation of the variograms are just done using the
#' spatial extent of the observations, not a possible temporal extent, such as
#' done by Skoien and Bloschl (2007). However, it is possible to make
#' predictions based on observations from different time steps, through the use
#' of the lag-vectors. These vectors describe a typical "delay" for each
#' observation and prediction location. This delay could for runoff related
#' variables be similar to travel time to each gauging location. For a certain
#' prediction location, earlier time steps would be picked for neighbours with
#' shorter travel time and later time steps for neighbours with slower travel
#' times.
#'
#' The lagExact parameter indicates whether to use a weighted average of two
#' time steps, or just the time step which is closest to the difference in lag
#' times.
#'
#' Universal kriging is invoked by giving \code{formulaString} a non-trivial
#' right hand side, e.g. \code{obs ~ elev} or \code{obs ~ x + y}, where the
#' trend variables are attribute columns of the observations and prediction
#' locations and/or the reserved coordinate names \code{x} and \code{y}. The
#' kriging system is then augmented with the trend basis functions instead of
#' only the unbiasedness constraint of ordinary kriging. The basis functions
#' are evaluated at the centroid of each support area, or block-averaged over
#' the discretisation points of \code{\link{rtopDisc}}, controlled by the
#' parameter \code{ukTrendSupport} (see \code{\link{getRtopParams}}). The
#' sample variogram is in this case computed from the residuals against an
#' OLS fit of the trend. For spatiotemporal interpolation the trend variables
#' have to be time-invariant columns of the \code{sp}-slot.
#'
#' The use of lag times should in theory increase the computation time, but
#' might, due to different computation methods, even speed up the computation
#' when the number of neighbours to be used (parameter nmax) is small compared
#' to the number of observations. If computation is slow, it can be useful to
#' test olags = rep(0, `dim(observations)[1]`) and similar for
#' predictionLocations.
#'
#' @param object object of class \code{rtop} or
#' \code{\link[sp]{SpatialPolygonsDataFrame}} or \code{\link[spacetime]{STSDF}}
#' @param varMatUpdate logical; if TRUE, also existing variance matrices will
#' be recomputed, if FALSE, only missing variance matrices will be computed,
#' see also \code{\link{varMat}}
#' @param predictionLocations \code{\link[sp]{SpatialPolygons}} or
#' \code{\link[sp]{SpatialPolygonsDataFrame}} or \code{\link[spacetime]{STSDF}}
#' with prediction locations. NULL if cross validation is to be performed.
#' @param varMatObs covariance matrix of observations, where diagonal must
#' consist of internal variance, typically generated from call to
#' \code{\link{varMat}}
#' @param varMatPredObs covariance matrix between observation locations and
#' prediction locations, typically generated from call to \code{\link{varMat}}
#' @param varMat list covariance matrices including the two above
#' @param params a set of parameters, used to modify the default parameters for
#' the \code{rtop} package, set in \code{\link{getRtopParams}}. Additionally,
#' it is possible overrule some of the parameters in \code{object$params} by
#' passing them as separate arguments.
#' @param formulaString formula that defines the dependent variable as a linear
#' model of independent variables, see e.g. \code{\link{createRtopObject}} for
#' more details.
#' @param sel array of prediction location numbers, if only a limited number of
#' locations are to be interpolated/crossvalidated
#' @param wret logical; if TRUE, return a matrix of weights instead of the
#' predictions, useful for batch processing of time series, see also details
#' @param dObs list of discretisation points for the observation areas from
#' \code{\link{rtopDisc}}, used for block-averaging of coordinate-based trend
#' basis functions when \code{ukTrendSupport = "block"}; computed internally
#' when needed and not supplied
#' @param dPred list of discretisation points for the prediction areas, see
#' \code{dObs}
#' @param trendObs matrix of trend basis functions evaluated at the
#' observation supports (one row per observation, including the intercept);
#' mainly for internal use, computed from \code{formulaString} when missing
#' @param trendPred matrix of trend basis functions evaluated at the
#' prediction supports, see \code{trendObs}
#' @param olags A vector describing the relative lag which should be applied
#' for the observation locations. See also details
#' @param plags A vector describing the relative lag which should be applied
#' for the predicitonLocations. See also details
#' @param lagExact logical; whether differences in lagtime should be computed
#' exactly or approximate
#' @param ... from \code{rtopKrige.rtop}, arguments to be passed to
#' \code{rtopKrige.default}. In \code{rtopKrige.default}, parameters for
#' modification of the object parameters or default parameters.  Of particular
#' interest are \code{cv}, a logical for doing cross-validation, \code{nmax},
#' and \code{maxdist} for maximum number of neighbours and maximum distance to
#' neighbours, respectively, and \code{wlim}, the limit for the absolute values
#' of the weights. It can also be useful to set \code{singularSolve} if some of
#' the areas are almost similar, see also details below.
#' @return If called with \code{\link[sp]{SpatialPolygonsDataFrame}}, the
#' function returns a \cr \code{\link[sp]{SpatialPolygonsDataFrame}} with
#' predictions, either at the locations defined in \cr
#' \code{predictionLocations}, or as leave-one-out cross-validation predicitons
#' at the same locations as in object if \code{cv = TRUE}
#'
#' If called with an rtop-object, the function returns the same object with the
#' predictions added to the object.
#' @author Jon Olav Skoien
#' @references Skoien J. O., R. Merz, and G. Bloschl. Top-kriging -
#' geostatistics on stream networks. Hydrology and Earth System Sciences,
#' 10:277-287, 2006.
#'
#' Skoien, J. O. and G. Bloschl. Spatio-Temporal Top-Kriging of Runoff Time
#' Series. Water Resources Research 43:W09419, 2007.
#'
#' Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione,
#' A., 2014. Rtop: An R package for interpolation of data with a variable
#' spatial support, with an example from river networks. Computers &
#' Geosciences, 67.
#' @keywords spatial
#' @examples
#'
#' \donttest{
#' # The following command will download  the complete example data set
#' # downloadRtopExampleData()
#' # observations$obs = observations$QSUMMER_OB/observations$AREASQKM
#'
#' rpath <- system.file("extdata",package="utop")
#' library(sf)
#' observations <- sf::st_read(rpath, "observations")
#' predictionLocations <- sf::st_read(rpath,"predictionLocations")
#'
#' # Setting some parameters; nclus > 1 will start a cluster with nclus
#' # workers for parallel processing
#' params <- list(gDist = TRUE, cloud = FALSE, nclus = 1, rresol = 25)
#'
#' # Create a column with the specific runoff:
#' observations$obs <- observations$QSUMMER_OB/observations$AREASQKM
#'
#' # Build an object
#' rtopObj <- createRtopObject(observations, predictionLocations,
#'                            params = params)
#'
#' # Fit a variogram (function also creates it)
#' rtopObj <- rtopFitVariogram(rtopObj)
#'
#' # Predicting at prediction locations
#' rtopObj <- rtopKrige(rtopObj)
#'
#' # Cross-validation
#' rtopObj <- rtopKrige(rtopObj,cv=TRUE)
#' cor(rtopObj$predictions$observed,rtopObj$predictions$var1.pred)
#' }
#'
#' @export
rtopKrige <- function(object, ...) UseMethod("rtopKrige")


#' Spatial simulation of data with spatial support
#'
#' rtopSim will conditionally or unconditionally simulate data with areal
#' support. This function should be seen as experimental, some issues are
#' described below.
#'
#' This function can do constrained or unconstrained simulation for areas. The
#' simplest way of calling the function is with an rtop-object that contains
#' the fitted variogram model and all the other necessary data (see
#' \code{\link{createRtopObject}} or \code{\link{utop-package}}).
#' \code{rtopSim} is the only function in \code{rtop} which does not need
#' observations. However, a variogram model is still necessary to perform
#' simulations.
#'
#' The arguments \code{beta} and \code{largeFirst} are only used for
#' unconditional simulations.
#'
#' The function is still in an experimental stage, and might change in the
#' future. There are some issues with the current implementation:
#' - Numerical issues can in some cases give negative estimation variances,
#'   which will result in an invalid distribution for the simulation. This will
#'   result in simulated NA values for these locations.
#' - The variability of simulated values for small areas (such as small
#'   headwater catchments) will be relatively high based on the statistical
#'   uncertainty. This could be overestimated compared to the uncertainty which
#'   is possible based on rainfall.
#'
#' @param object object of class \code{rtop} or
#' \code{\link[sp]{SpatialPolygonsDataFrame}} or \code{sf}
#' (\code{\link[sf]{st_sf}}) or \code{NULL}
#' @param varMatUpdate logical; if TRUE, also existing variance matrices will
#' be recomputed, if FALSE, only missing variance matrices will be computed,
#' see also \code{\link{varMat}}
#' @param beta The expected mean of the data, for unconditional simulations
#' @param largeFirst Although the simulation method follows a random path
#' around the predictionLocations, simulating the largest area first will
#' assure that the true mean of the simulated values will be closer to beta
#' @param replace logical; if observation locations are also present as
#' predictions, should they be replaced?  This is particularly when doing
#' conditional simulations for a set of observations with uncertainty.
#' @param params a set of parameters, used to modify the standard parameters
#' for the \code{rtop} package, set in \code{\link{getRtopParams}}. The
#' argument params can also be used for the other methods, through the
#' ...-argument.
#' @param dump file name for saving the updated object, after adding variance
#' matrices. Useful if there are problems with the simulation, particularly if
#' it for some reason crashes.
#' @param debug.level logical that controls some output, will override the
#' object parameters
#' @param predictionLocations \code{\link[sp]{SpatialPolygons}} or
#' \code{\link[sp]{SpatialPolygonsDataFrame}} or \code{sf}-object with
#' locations for simulations.
#' @param varMatObs covariance matrix of possible observations, where diagonal
#' must consist of internal variance, typically generated from call to
#' \code{\link{varMat}}
#' @param varMatPredObs covariance matrix between possible observation
#' locations and simulation locations, typically generated from call to
#' \code{\link{varMat}}
#' @param varMatPred covariance matrix between simulation locations, typically
#' generated from a call to \code{\link{varMat}}
#' @param variogramModel a variogram model of type
#' \code{\link{rtopVariogramModel}}
#' @param ... possible modification of the object parameters or default
#' parameters.
#' @return If called with \code{\link[sp]{SpatialPolygons}} or \code{sf} as
#' predictionLocations and either \cr
#' \code{\link[sp]{SpatialPolygonsDataFrame}}, \code{sf} or \code{NULL} for
#' observations, the function returns a\cr
#' \code{\link[sp]{SpatialPolygonsDataFrame}} or \code{sf} with simulations at
#' the locations defined in \cr \code{predictionLocations}
#'
#' If called with an rtop-object, the function returns the same object with the
#' simulations added to the object.
#' @author Jon Olav Skoien
#' @references Skoien J. O., R. Merz, and G. Bloschl. Top-kriging -
#' geostatistics on stream networks. Hydrology and Earth System Sciences,
#' 10:277-287, 2006.
#'
#' Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione,
#' A., 2014. Rtop: An R package for interpolation of data with a variable
#' spatial support, with an example from river networks. Computers &
#' Geosciences, 67.
#' @keywords spatial
#' @examples
#'
#' \donttest{
#' # The following command will download  the complete example data set
#' # downloadRtopExampleData()
#'
#' rpath <- system.file("extdata",package="utop")
#' library(sf)
#' observations <- sf::st_read(rpath, "observations")
#' predictionLocations <- sf::st_read(rpath,"predictionLocations")
#'
#' # Setting some parameters; nclus > 1 will start a cluster with nclus
#' # workers for parallel processing
#' params <- list(gDist = TRUE, cloud = FALSE, nclus = 1, rresol = 25)
#'
#' # Create a column with the specific runoff:
#' observations$obs <- observations$QSUMMER_OB/observations$AREASQKM
#'
#' # Build an object
#' rtopObj <- createRtopObject(observations, predictionLocations,
#'                            params = params, formulaString = "obs~1")
#'
#' # Fit a variogram (function also creates it)
#' rtopObj <- rtopFitVariogram(rtopObj)
#'
#' # Conditional simulations for two new locations
#' rtopObj10 <- rtopSim(rtopObj, nsim = 5)
#' rtopObj11 <- rtopObj
#'
#' # Unconditional simulation at the observation locations
#' # (These are moved to the predictionLocations)
#' rtopObj11$predictionLocations <- rtopObj11$observations
#' rtopObj11$observations <- NULL
#' # Setting varMatUpdate to TRUE, to make sure that covariance
#' # matrices are recomputed
#' rtopObj12 <- rtopSim(rtopObj11, nsim = 10, beta = 0.01,
#'                     varMatUpdate = TRUE)
#'
#' sp::summary(data.frame(rtopObj10$simulations))
#' sp::summary(data.frame(rtopObj12$simulations))
#'
#' }
#'
#' @export
rtopSim <- function(object, ...) UseMethod("rtopSim")
#' @rdname rtopVariogramModel
#' @export
updateRtopVariogram <- function(object, ...) UseMethod("updateRtopVariogram")
