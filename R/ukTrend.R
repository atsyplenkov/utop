# Universal-kriging trend helpers.
#
# The trend (drift) is defined by the RHS of formulaString. Basis functions
# can be attributes of the areas and/or the reserved coordinate names
# "x"/"y", which are filled from the support geometry. With areal support,
# basis functions that vary within an area (coordinate-based terms) are
# either evaluated at the centroid or block-averaged over the rtopDisc()
# discretisation points, controlled by params$ukTrendSupport.

#' @noRd
ukTrendTerms <- function(formulaString) {
  if (!inherits(formulaString, "formula")) {
    formulaString <- as.formula(formulaString)
  }
  delete.response(terms(formulaString))
}

#' @noRd
hasUkTrend <- function(formulaString) {
  length(attr(ukTrendTerms(formulaString), "term.labels")) > 0
}

#' @noRd
ukNLoc <- function(locations) {
  if (inherits(locations, "stars")) {
    utop_stars_nspace(locations)
  } else {
    nrow(utop_as_sf(locations))
  }
}

#' @noRd
#' @noRd
ukLocData <- function(locations) {
  sf::st_drop_geometry(utop_as_sf(locations))
}

#' @noRd
#' @noRd
ukCentroids <- function(locations) {
  utop_centroid_coordinates(locations)
}

#' @noRd
#' @noRd
ukDiscCoordinates <- function(discPoints) {
  utop_point_coordinates(discPoints)
}

#' @noRd
ukModelMatrix <- function(tt, data) {
  mf <- model.frame(tt, data, na.action = na.fail)
  isnum <- vapply(
    mf,
    FUN = function(v) is.numeric(v) || is.logical(v),
    FUN.VALUE = logical(1)
  )
  if (!all(isnum)) {
    stop(paste(
      "Universal kriging trend variables must be numeric:",
      paste(names(mf)[!isnum], collapse = ", ")
    ))
  }
  model.matrix(tt, mf)
}

#' Evaluate the trend basis functions of the RHS of formulaString for a set
#' of support areas, either at centroids or block-averaged over the
#' discretisation points from rtopDisc()
#' @noRd
ukTrendMatrix <- function(
  formulaString,
  locations,
  params = NULL,
  discPoints = NULL
) {
  params <- coerce_utop_params(params)
  tt <- ukTrendTerms(formulaString)
  nloc <- ukNLoc(locations)
  if (length(attr(tt, "term.labels")) == 0) {
    if (!attr(tt, "intercept")) {
      stop("formulaString without intercept needs at least one trend variable")
    }
    return(matrix(
      1,
      nrow = nloc,
      ncol = 1,
      dimnames = list(NULL, "(Intercept)")
    ))
  }
  support <- params@uk_trend_support
  if (!support %in% c("centroid", "block")) {
    stop(paste(
      "Unknown ukTrendSupport:",
      support,
      "- must be \"centroid\" or \"block\""
    ))
  }
  df <- ukLocData(locations)
  vars <- all.vars(tt)
  coordVars <- intersect(setdiff(vars, names(df)), c("x", "y"))
  missVars <- setdiff(vars, c(names(df), coordVars))
  if (length(missVars) > 0) {
    stop(paste(
      "Universal kriging trend variables not found in data:",
      paste(missVars, collapse = ", "),
      "- only attribute columns and the coordinate names x/y can be used"
    ))
  }

  if (support == "block" && length(coordVars) > 0) {
    # Coordinate-based basis functions vary within each area; block-average
    # them over the discretisation points. Pure attribute terms are constant
    # within an area and are handled by the centroid branch below.
    if (is.null(discPoints)) {
      if (!inherits(locations, "sf") && !inherits(locations, "stars")) {
        stop(paste(
          "ukTrendSupport = \"block\" with coordinate-based trend terms",
          "requires polygon supports or precomputed discretisation points"
        ))
      }
      loc_sf <- utop_as_sf(locations)
      if (!all(sf::st_dimension(loc_sf) == 2)) {
        stop(paste(
          "ukTrendSupport = \"block\" with coordinate-based trend terms",
          "requires polygon supports or precomputed discretisation points"
        ))
      }
      discPoints <- utop_disc(locations, params = params)
    }
    if (length(discPoints) != nloc) {
      stop("Number of discretisation elements does not match number of areas")
    }
    Fmat <- NULL
    for (iloc in seq_len(nloc)) {
      pts <- ukDiscCoordinates(discPoints[[iloc]])
      ldf <- df[rep(iloc, dim(pts)[1]), , drop = FALSE]
      if ("x" %in% coordVars) {
        ldf$x <- pts[, 1]
      }
      if ("y" %in% coordVars) {
        ldf$y <- pts[, 2]
      }
      mm <- ukModelMatrix(tt, ldf)
      if (is.null(Fmat)) {
        Fmat <- matrix(
          0,
          nrow = nloc,
          ncol = dim(mm)[2],
          dimnames = list(NULL, colnames(mm))
        )
      }
      Fmat[iloc, ] <- colMeans(mm)
    }
    Fmat
  } else {
    if (length(coordVars) > 0) {
      cors <- ukCentroids(locations)
      if ("x" %in% coordVars) {
        df$x <- cors[, 1]
      }
      if ("y" %in% coordVars) {
        df$y <- cors[, 2]
      }
    }
    Fmat <- ukModelMatrix(tt, df)
    dimnames(Fmat) <- list(NULL, colnames(Fmat))
    Fmat
  }
}

#' OLS residuals of the dependent variable against the trend basis functions,
#' used for the residual sample variogram
#' @noRd
ukResiduals <- function(
  formulaString,
  locations = NULL,
  params = list(),
  discPoints = NULL,
  depValues = NULL,
  trendMatrix = NULL
) {
  if (is.null(trendMatrix)) {
    trendMatrix <- ukTrendMatrix(formulaString, locations, params, discPoints)
  }
  if (is.null(depValues)) {
    if (!inherits(formulaString, "formula")) {
      formulaString <- as.formula(formulaString)
    }
    depValues <- ukLocData(locations)[[as.character(formulaString[[2]])]]
  }
  ok <- which(!is.na(depValues))
  beta <- qr.coef(qr(trendMatrix[ok, , drop = FALSE]), depValues[ok])
  beta[is.na(beta)] <- 0
  as.vector(depValues - trendMatrix %*% beta)
}
