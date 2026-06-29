# For compiling the Fortran file:
# R CMD SHLIB vred.f

#' Create an object for interpolation within the utop package
#'
#' This is a help function for creating an object (see
#' \code{\link{utop-package}} to be used for interpolation within the utop
#' package
#'
#'
#' @param observations \code{\link[sf]{sf}} polygons, or a vector data cube
#'   with \code{stars} for spatiotemporal observations
#' @param predictionLocations \code{\link[sf]{sf}} polygons, or a vector data
#'   cube with \code{stars} for spatiotemporal prediction
#'   locations
#' @param formulaString formula that defines the dependent variable as a linear
#' model of independent variables; suppose the dependent variable has name
#' \code{z}, for ordinary and simple kriging use the formula \code{z~1}; for
#' universal kriging, suppose \code{z} is linearly dependent on \code{x} and
#' \code{y}, use the formula \code{z~x+y}. The formulaString defaults to
#' \code{"value~1"} if \code{value} is a part of the data set.  If not, the
#' first column of the data set is used. The trend variables on the RHS can
#' be attribute columns of the observations and prediction locations and/or
#' the reserved coordinate names \code{x} and \code{y}. The trend basis
#' functions are evaluated either at the centroids of the areas or
#' block-averaged over the discretisation points from \code{\link{rtopDisc}},
#' controlled by the parameter \code{ukTrendSupport}, see
#' \code{\link{getRtopParams}}. For stars vector cubes, RHS covariates must
#' be time-invariant for each spatial support.
#' @param params parameters to modify the default parameters of the
#' utop-package, set internally in this function by a call to
#' \code{\link{getRtopParams}}
#' @param overlapObs matrix with observations that overlap each other
#' @param overlapPredObs matrix with \code{observations} and
#' \code{predictionLocations} that overlap each other
#' @param ... Extra parameters to \code{\link{getRtopParams}} and possibility
#' to pass depreceted arguments
#' @return An object of class \code{rtop} with observations, prediction
#' locations, parameters and possible other elements useful for interpolation
#' in the utop-package. Most other externally visible functions in the
#' package will be able to work with this object, and add the results as a new
#' element.
#' @author Jon Olav Skoien
#' @seealso \code{\link{getRtopParams}}
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
#' rtopObj <- createRtopObject(observations, predictionLocations,
#'                            params = params)
#' }
#'
#'
#' @export
createRtopObject <- function(
  observations,
  predictionLocations,
  formulaString,
  params = list(),
  overlapObs,
  overlapPredObs,
  ...
) {
  dots <- list(...)
  if (inherits(observations, "rtop")) {
    # Updating object with parameters
    object <- observations
    object$params <- getRtopParams(
      params,
      observations = object$observations,
      formulaString = if (missing(formulaString)) {
        object$formulaString
      } else {
        formulaString
      },
      ...
    )
    return(object)
  }
  object <- list()

  if (missing(observations)) {
    stop("Observations are missing")
  }
  if (utop_is_legacy_sp(observations)) {
    utop_stop_legacy_sp(observations)
  }
  observations <- utop_add_area(observations)

  object$observations <- observations

  if (!missing(predictionLocations)) {
    if (utop_is_legacy_sp(predictionLocations)) {
      utop_stop_legacy_sp(predictionLocations)
    }
    predictionLocations <- utop_add_area(predictionLocations)

    obs_crs <- sf::st_crs(utop_as_sf(observations))
    pred_crs <- sf::st_crs(utop_as_sf(predictionLocations))
    if (!isTRUE(all.equal(is.na(obs_crs), is.na(pred_crs)))) {
      stop("only one of observations and predictionLocations have projection")
    }
    if (!is.na(obs_crs) && obs_crs != pred_crs) {
      stop(paste(
        "observations and predictionLocations have different projections:",
        obs_crs,
        pred_crs
      ))
    }

    object$predictionLocations <- predictionLocations
  }
  if (missing(formulaString)) {
    formulaString <- utop_default_formula(observations)
    warning(paste("formulaString missing, using", formulaString))
  }
  if (!inherits(formulaString, "formula")) {
    formulaString <- as.formula(formulaString)
  }
  object$formulaString <- formulaString
  #  depVar = formulaString[[2]] else depVar = "obs"
  object$params <- getRtopParams(
    newPar = params,
    formulaString = formulaString,
    observations = observations
  )
  if (length(dots) > 0) {
    object <- modifyList(object, dots)
  }
  if (object$params$nugget) {
    if (!missing(overlapObs) && !is.null(overlapObs)) {
      object$overlapObs <- overlapObs
    } else {
      object$overlapObs <- findOverlap(
        observations,
        debug.level = object$params$debug.level
      )
    }
    if (!missing(overlapPredObs) && !is.null(overlapPredObs)) {
      object$overlapPredObs <- overlapPredObs
    } else if (!missing(predictionLocations)) {
      object$overlapPredObs <- findOverlap(
        observations,
        predictionLocations,
        debug.level = object$params$debug.level
      )
    }
  }
  class(object) <- "rtop"
  object
}


#' Setting parameters for the intamap package
#'
#' This function sets a range of the parameters for the intamap package, to be
#' included in the object described in \code{\link{utop-package}}
#'
#'
#' @param params An existing set of parameters for the interpolation process,
#' of class \cr \code{intamapParams} or a list of parameters for modification
#' of the default parameters
#' @param newPar A \code{list} of parameters for updating \code{params} or for
#' modification of the default parameters.  Possible parameters with their
#' defaults are given below
#' @param observations \code{\link[sf]{sf}} polygons, used for setting some
#' of the default parameters
#' @param formulaString formula that defines the dependent variable as a linear
#' model of independent variables, see e.g. \code{\link{createRtopObject}} for
#' more details.
#' @param ... Individual parameters for updating \code{params} or for
#' modification of the default parameters.  Possible parameters with their
#' defaults are given below
#'
#' \describe{ \item{model = "Ex1"}{ - variogram model type. Currently the
#' following models are implemented: \describe{ \item{Exp}{ - Exponential
#' model} \item{Ex1}{ - Multiplication of a modified exponential and fractal
#' model, the same model as used in Skoien et al(2006).} \item{Gau}{ - Gaussian
#' model} \item{Ga1}{ - Multiplication of gaussian and fractal model}
#' \item{Sph}{ - Spherical model} \item{Sp1}{ - Multiplication of spherical and
#' fractal model} \item{Fra}{ - Fractal model} }} \item{parInit}{ - the initial
#' parameters and the limits of the variogram model to be fitted, given as a
#' matrix with three columns, where the first column is the lower limit, the
#' second column is the upper limit and the third column are starting values.}
#' \item{nugget = FALSE}{ - logical; if TRUE, nugget effect should be
#' estimated} \item{unc = TRUE}{ - logical; if TRUE the variance of
#' observations are in column \code{unc}} \item{rresol = 100}{ - minimum number
#' of discretization points in each area} \item{hresol = 5}{ - number of
#' discretization points in one direction for elements in binned variograms}
#' \item{cloud = FALSE}{ - logical; if TRUE use the cloud variogram for
#' variogram fitting} \item{amul = 1}{ - defines the number of areal bins
#' within one order of magnitude. Numbers between 1 and 3 are possible, as this
#' parameter refers to the \code{axp} parameter of
#' \code{\link[graphics]{axTicks}}.} \item{dmul = 3}{ - defines the number of
#' distance bins within one order of magnitude. Numbers between 1 and 3 are
#' possible, as this parameter refers to the \code{axp} parameter of
#' \code{\link[graphics]{axTicks}}.} \item{fit.method = 9}{ - defines the type
#' of Least Square method for fitting of variogram.  The methods 1-7 correspond
#' to the similar methods in \code{\link[gstat]{fit.variogram}}.  \describe{
#' \item{1}{ - weighted least squares with number of pairs per bin: \cr err = n *
#' (yobs-ymod)^2} \item{2}{ - weighted least squares difference according to
#' Cressie (1985): \cr err2 = abs(yobs/ymod-1)} \item{6}{ - ordinary least
#' squares difference: err = (yobs-ymod)^2} \item{7}{ - similar to default of
#' gstat, where higher weights are given to shorter distances err = n/h^2 *
#' (yobs-mod)^2} \item{8}{ - Opposite of weighted least squares difference
#' according to Cressie (1985): err3=abs(ymod/yobs-1)} \item{9}{ - neutral
#' WLS-method - err = min(err2,err3)} } } \item{gDistEst = FALSE}{ - use
#' geostatistical distance when fitting variograms} \item{gDistPred = FALSE}{ -
#' use geostatistical distance for semivariogram matrices} \item{gDist}{ -
#' parameter to set jointly \code{gDistEst = gDistPred = gDist}} \item{nmax =
#' 10}{for local kriging: the number of nearest observations that should be
#' used for a kriging prediction or simulation, where nearest is defined in
#' terms of the space of the spatial locations.  By default, 10 observations
#' are used.} \item{maxdist = Inf}{ - for local kriging: only observations
#' within a distance of \code{maxdist} from the prediction location are used
#' for prediction or simulation; if combined with nmax, both criteria apply }
#' \item{hstype = "regular"}{ - sampling type for binned variograms}
#' \item{rstype = "rtop"}{ - sampling type for the elements, see also
#' \code{\link{rtopDisc}}} \item{nclus = 1}{- number of CPUs to use if parallel
#' processing is wanted; nclus = 1 means no parallelization} \item{cnAreas =
#' 100}{- limit whether parallel processing should be applied; the minimum
#' number of areas in \code{\link{varMat}}, and also controlling when to use
#' parallel processing in \code{\link{rtopDisc}}, when \cr
#' \code{nAreas*params$rresol/100 > cnAreas} } \item{clusType = NULL}{- the
#' cluster type to be started for parallel processing; uses the default type of
#' the system when clusType = NULL} \item{outfile = NULL}{file where output can
#' be printed during parallel execution} \item{varClean = FALSE}{logical; if
#' TRUE it will remove highly correlated areas from the covariance matrix
#' during simulation } \item{wlim = 1.5}{ - an upper limit for the norm of the
#' weights in kriging, see \code{\link{rtopKrige}}} \item{wlimMethod =
#' "all"}{which method to use for reducing the norm of the weights if
#' necessary. Either "all", which modifies all weights equally or "neg" which
#' reduces negative weights and large weights more than the smallest weights }
#' \item{singularSolve}{ - logical; When TRUE, the kriging function will
#' attempt to solve singular kriging matrices by removing catchments that have
#' the same correlations. This will usually happen when two catchments are
#' almost overlapping, and they are discretized with the same points. See also
#' \code{\link{rtopKrige}}.} \item{ukTrendSupport = "centroid"}{ - how the
#' universal kriging trend basis functions (the RHS of \code{formulaString})
#' are evaluated for each support area. \code{"centroid"} evaluates them at
#' the area centroid, \code{"block"} averages them over the discretisation
#' points of the area from \code{\link{rtopDisc}}. The difference only
#' matters for basis functions that vary within an area, i.e., terms using
#' the reserved coordinate names \code{x} and \code{y}; attribute covariates
#' are constant within an area.} \item{cv = FALSE}{ - logical; for cross-validation
#' of observations} \item{debug.level = 1}{ - used in some functions for giving
#' additional output. See individual functions for more information.}
#' \item{partialOverlap = FALSE}{whether to work with partially overlapping
#' areas} \item{olim = 1e-4}{smallest overlapping area to be used for partial
#' overlap, relative to the smallest of the areas} \item{nclus = 1}{option to
#' use parallel processing, nclus > 1 defines the number of workers to be
#' started} \item{clusType = NA}{which cluster type to start if nclus > 1; the
#' default is used if nclusType = NA } \item{cnAreas = 200}{The minimum number
#' of observations or observations plus predictions allowing parallelization in
#' the creation of the covariance matrix} \item{cDlim = 1e6}{The minimum number
#' of discretization points for allowing parallelization in the discretization
#' process} \item{observations}{ - used for initial values of parameters if
#' supplied} \item{formulaString}{ - used for initial values of parameters if
#' supplied} }
#' @return A list of the parameters with class \code{rtopParams} to be included
#' in the \code{object} described in \link{utop-package}
#' @note This function will mainly be called by \code{\link{createRtopObject}},
#' but can also be called by the user to create a parameter set or update an
#' existing parameter set. If none of the arguments is a list of class
#' \code{rtopParams}, the function will assume that the argument(s) are
#' modifications to the default set of parameters. The function can also be
#' called by other functions in the utop-package if the users chooses not to
#' work with an object of class \code{rtop}.
#'
#' If the function is called with two lists of parameters (but the first one is
#' not of class \code{rtopParams}) they are both seen as modifications to the
#' default parameter set. If they share some parameters, the parameter values
#' from the second list will be applied.
#'
#' Parallel processing has been included for some of the functions. The default
#' is no parallel procesing, and the package also attempts to decide whether it
#' is sensible to start a set of clusters and distribute jobs to them based on
#' the size of the job. The default limit might not be the best for every
#' system.
#' @author Jon Olav Skoien
#' @seealso \code{\link{createRtopObject}} and \code{\link{utop-package}}
#' @references Cressie, N. 1985. Fitting variogram models by weighted least
#' squares. Mathematical Geology, 17 (5), 563-586
#'
#' Skoien J. O., R. Merz, and G. Bloschl. Top-kriging - geostatistics on stream
#' networks. Hydrology and Earth System Sciences, 10:277-287, 2006
#'
#' Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione,
#' A., 2014. Rtop: An R package for interpolation of data with a variable
#' spatial support, with an example from river networks. Computers &
#' Geosciences, 67.
#' @keywords spatial
#' @examples
#'
#' # Create a new set of intamapParameters, with default parameters:
#' params <- getRtopParams()
#' # Make modifications to the default list of parameters
#' params <- getRtopParams(newPar = list(gDist = TRUE, nugget = FALSE))
#' # Make modifications to an existing list of parameters
#' params <- getRtopParams(params = params, newPar = list(gDist = TRUE,
#'          nugget = FALSE))
#'
#' @export
getRtopParams <- function(
  params = list(),
  newPar = list(),
  observations,
  formulaString,
  ...
) {
  oldPar <- params
  dots <- list(...)
  if (inherits(oldPar, "intamapParams") || inherits(newPar, "intamapParams")) {
    intPar <- TRUE
  } else {
    intPar <- FALSE
  }
  oClass <- class(oldPar)
  nClass <- class(newPar)
  if (inherits(oldPar, "rtopParams")) {
    params <- oldPar
    oldPar <- list()
  } else if (inherits(newPar, "rtopParams")) {
    params <- newPar
    newPar <- list()
  } else {
    params <- getRtopDefaultParams(...)
  }
  if (
    length(grep("geoDist", names(oldPar))) > 0 ||
      length(grep("geoDist", names(newPar))) > 0 ||
      length(grep("geoDist", names(dots))) > 0
  ) {
    stop("geoDist is not used anymore, please use gDist")
  }

  params <- modifyList(params, oldPar)
  params <- modifyList(params, newPar)
  gDist <- ifelse(
    "gDist" %in% names(dots),
    dots$gDist,
    ifelse(
      "gDist" %in% names(newPar),
      newPar$gDist,
      ifelse("gDist" %in% names(oldPar), oldPar$gDist, FALSE)
    )
  )
  if (gDist) {
    params$gDistEst <- TRUE
    params$gDistPred <- TRUE
  }

  if (!missing(observations) && !("parInit" %in% names(params))) {
    if (missing(formulaString)) {
      formulaString <- utop_default_formula(
        observations,
        caller = "getRtopParams: formulaString"
      )
      warning(paste(
        "getRtopParams: formulaString missing, using",
        formulaString
      ))
    }
    params$parInit <- findParInit(formulaString, observations, params$model)
  } else if (!("parInit" %in% names(params))) {
    params$parInit <- findParInitDefault(params$model)
  }
  params <- modifyList(params, dots)

  if (intPar) {
    class(params) <- c("rtopParams", "intamapParams")
  } else {
    class(params) <- "rtopParams"
  }
  params
}


#' @noRd
getRtopDefaultParams <- function(
  parInit,
  model = "Ex1",
  nugget = FALSE,
  unc = TRUE,
  rresol = 100, # Resolution real areas
  hresol = 5, # Resolution in x-direction rectangles
  #   logtrans = FALSE, # Logtransform data
  cloud = FALSE, # work with cloud variogram
  #   cutoff,        # cutoff distance in variogram - better to set in ... in call to function,
  amul = 2, # amul - defines the number of areal bins within one order of magnitude
  dmul = 3, # dmul - defines the number of distance bins within one order of magnitude
  fit.method = 9, # ils - Defines the type of Least Square method for fitting of variogram
  #       1 - least squares difference  - err = yobs-ymod
  #       2 - Weighted least squares difference according to Cressie (1985) - err2=n(yobs/ymod-1)^2
  #       6 - No weights
  #       7 - gstat fitting (Nj/hj^2)
  #       8 - opposite of weighted least squares difference according to Cressie (1985) - err2=n*(ymod/yobs-1)^2
  #       9 - Neutreal WLS-method - err = min(err2,err3)
  gDistEst = FALSE, # use ghosh distance
  gDistPred = FALSE,
  varClean = FALSE,
  maxdist = Inf,
  nmax = 10,
  hstype = "regular", # Sampling type for hypothetical areas
  # Sampling type for real areas
  rstype = "rtop",
  nclus = 1,
  cnAreas = 100,
  clusType = NULL,
  outfile = NULL,
  partialOverlap = FALSE,
  wlim = 1.5,
  wlimMethod = "all",
  singularSolve = FALSE,
  ukTrendSupport = "centroid",
  cv = FALSE,
  debug.level = if (interactive()) 1 else 0,
  observations,
  formulaString
) {
  list(
    model = model,
    nugget = nugget,
    unc = unc,
    rresol = rresol,
    hresol = hresol,
    rstype = rstype,
    hstype = hstype,
    #     logtrans = logtrans,
    cloud = cloud,
    #     cutoff = cutoff,
    amul = amul,
    dmul = dmul,
    fit.method = fit.method,
    gDistEst = gDistEst,
    gDistPred = gDistPred,
    varClean = varClean,
    maxdist = maxdist,
    nmax = nmax,
    nclus = nclus,
    cnAreas = cnAreas,
    clusType = clusType,
    outfile = outfile,
    partialOverlap = partialOverlap,
    wlim = wlim,
    wlimMethod = wlimMethod,
    singularSolve = singularSolve,
    ukTrendSupport = ukTrendSupport,
    cv = cv,
    debug.level = debug.level
  )
}

###########################################
#' @noRd
findParInitDefault <- function(model) {
  #  parameters are: sill, range, nugget, fractal, weibull par
  parInit <- data.frame(
    parl = c(1e-06, 1e-02, 1.0e-01, 1e-5, 1e-01),
    paru = c(5.0e+02, 1.0e7, 1.0e+07, 1.5, 1.7)
  )
  parInit$par0 <- 10**(0.5 * (log10(parInit$paru) + log10(parInit$parl)))

  if (model %in% c("Exp", "Sph", "Gau")) {
    parInit <- parInit[1:3, ]
  } else if (model == "Sp1") {
    parInit <- parInit[1:4, ]
  } else if (model == "Ex1") {
    parInit <- parInit
  } else if (model == "Fra") {
    parInit[2, ] <- c(1e-6, 2, 0.01)
    parInit <- parInit[1:3, ]
  } else {
    stop(paste("model", model, "not implemented"))
  }
  parInit
}

#########################################
#' @noRd
findParInit <- function(formulaString, observations, model) {
  observations <- utop_add_area(observations)
  if (inherits(observations, "stars")) {
    ntime <- utop_stars_ntime(observations)
    if (ntime > 20) {
      observations <- utop_stars_slice_time(
        observations,
        sample(seq_len(ntime), 20)
      )
    }
    vario <- rtopVariogram(observations, formulaString = formulaString)
    aObs <- utop_area(observations)
  } else {
    if (hasUkTrend(formulaString)) {
      observations$ukResidual <- ukResiduals(formulaString, observations)
      formulaString <- ukResidual ~ 1
    }
    vario <- gstat::variogram(formulaString, observations)
    aObs <- observations$area
  }
  parInit <- data.frame(parl = c(1:5), paru = 1, par0 = 1)
  parInit[1, 1] <- min(vario$gamma) / 10
  parInit[1, 2] <- max(vario$gamma) * 500
  parInit[2, 1] <- sqrt(min(aObs)) / 4
  parInit[2, 2] <- max(vario$dist) * 10
  minla <- min(aObs)
  maxla <- (max(aObs)^1.5) * max(vario$gamma)
  parInit[3, 1] <- min(vario$gamma) * minla / 100
  parInit[3, 2] <- max(vario$gamma) * maxla
  parInit[4, 1] <- 1e-5
  parInit[4, 2] <- 1.5
  parInit[5, 1] <- 0.1
  parInit[5, 2] <- 1.7
  if (model == "Ex1") {
    parInit[4, 2] <- 1
    parInit[5, 2] <- 1
  }

  parInit[, 3] <- sqrt(parInit[, 1] * parInit[, 2])
  if (model %in% c("Exp", "Sph", "Gau")) {
    parInit <- parInit[1:3, ]
  } else if (model == "Sp1") {
    parInit <- parInit[1:4, ]
  } else if (model == "Ex1") {
    parInit <- parInit
  } else if (model == "Fra") {
    parInit[2, ] <- c(1e-6, 2, 0.01)
    parInit <- parInit[1:3, ]
  } else {
    stop(paste("model", model, "not implemented"))
  }
  parInit
}
