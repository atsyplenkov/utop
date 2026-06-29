#' @export
#' @rdname rtopKrige
rtopKrige.rtop <- function(object, varMatUpdate = FALSE, params = list(), ...) {
  params <- getRtopParams(object$params, newPar = params, ...)
  if (!is.null(params$nsim) && params$nsim > 0) {
    return(rtopSim(object, varMatUpdate, params = params))
  }
  observations <- object$observations

  predictionLocations <- object$predictionLocations
  if (
    !all(c("varMatObs", "varMatPredObs") %in% names(object)) || varMatUpdate
  ) {
    object <- varMat(object, varMatUpdate, params = params, ...)
  }

  varMatObs <- object$varMatObs
  varMatPredObs <- object$varMatPredObs

  krigeRes <- rtopKrige(
    object = observations,
    predictionLocations = predictionLocations,
    varMatObs = varMatObs,
    varMatPredObs = varMatPredObs,
    params = params,
    formulaString = object$formulaString,
    dObs = object$dObs,
    dPred = object$dPred,
    ...
  )
  object$predictions <- krigeRes$predictions
  if ("cvInfo" %in% names(krigeRes)) {
    object$cvInfo <- krigeRes$cvInfo
  }
  if ("weight" %in% names(krigeRes)) {
    object$weight <- krigeRes$weight
  }
  if ("removed" %in% names(krigeRes)) {
    object$removed <- krigeRes$removed
  }
  object
}


#' @export
#' @noRd
rtopKrige.SpatialPolygonsDataFrame <- function(object, ...) {
  utop_stop_legacy_sp(object)
}


#' @export
#' @rdname rtopKrige
rtopKrige.default <- function(
  object,
  predictionLocations = NULL,
  varMatObs,
  varMatPredObs,
  varMat,
  params = list(),
  formulaString,
  sel,
  wret = FALSE,
  dObs = NULL,
  dPred = NULL,
  trendObs = NULL,
  trendPred = NULL,
  ...
) {
  params <- getRtopParams(params, ...)
  if (!is.null(params$nsim) && params$nsim > 0) {
    return(rtopSim(
      object,
      predictionLocations,
      varMatObs,
      varMatPredObs,
      varMat,
      params,
      formulaString,
      ...
    ))
  }
  #
  cv <- params$cv
  #  else object$params$cv = params$cv = cv
  nmax <- params$nmax
  wlim <- params$wlim
  wlimMethod <- params$wlimMethod
  maxdist <- params$maxdist
  debug.level <- params$debug.level
  lambda <- params$lambda
  if ("singularSolve" %in% names(params)) {
    singularSolve <- params$singularSolve
  } else {
    singularSolve <- FALSE
  }
  if (!is.null(lambda)) {
    BLUE <- TRUE
  } else {
    BLUE <- FALSE
  }
  if (!missing(varMat)) {
    if (missing(varMatObs)) {
      varMatObs <- varMat$varMatObs
    }
    if ((missing(varMatPredObs) | is.null(varMatPredObs)) && !cv) {
      varMatPredObs <- varMat$varMatPredObs
    }
  }
  depVar <- as.character(formulaString[[2]])
  observations <- object
  obs0 <- observations[[depVar]]
  nobs <- dim(object)[1]
  obscors <- utop_centroid_coordinates(observations)
  if (cv) {
    newcors <- obscors
  } else {
    newcors <- utop_centroid_coordinates(predictionLocations)
  }

  npred <- ifelse(cv, nobs, ifelse(!missing(sel), length(sel), dim(newcors)[1]))
  if (missing(sel)) {
    sel <- c(1:npred)
  }
  if (params$unc && "unc" %in% names(observations)) {
    unc0 <- observations$unc
  } else {
    unc0 <- array(0, nobs)
  }
  #
  # Universal kriging trend basis functions from the RHS of formulaString,
  # evaluated at centroids or block-averaged over the rtopDisc()
  # discretisation, according to params$ukTrendSupport.
  if (is.null(trendObs)) {
    trendObs <- ukTrendMatrix(formulaString, observations, params, dObs)
  }
  if (is.null(trendPred)) {
    if (cv) {
      trendPred <- trendObs
    } else {
      trendPred <- ukTrendMatrix(
        formulaString,
        predictionLocations,
        params,
        dPred
      )
    }
  } else if (!is.matrix(trendPred)) {
    trendPred <- matrix(trendPred, nrow = 1)
  }
  ptrend <- dim(trendObs)[2]
  if (dim(trendPred)[2] != ptrend) {
    stop(paste(
      "Different number of trend basis functions for observations",
      "and prediction locations"
    ))
  }
  #
  mdist <- sqrt(bbArea(utop_bbox(observations)))
  if (nobs < nmax && mdist < maxdist && !cv) {
    varMat <- rbind(varMatObs, t(trendObs))
    varMat <- cbind(varMat, rbind(trendObs, matrix(0, ptrend, ptrend)))
    diag(varMat)[1:nobs] <- unc0
    varInv <- solve(varMat)
    singMat <- TRUE
  } else {
    singMat <- FALSE
  }
  #
  if (cv) {
    predictionLocations <- observations
    varMatPredObs <- varMatObs
  }
  #
  removed <- NULL
  predictions <- data.frame(
    var1.pred = rep(0, npred),
    var1.var = 0,
    sumWeights = 0
  )
  if (cv) {
    predictions <- cbind(
      predictions,
      observed = observations[[depVar]],
      residual = 0,
      zscore = 0
    )
    cvInfo0 <- data.frame(
      inew = c(0),
      obs = c(1),
      pred = c(0),
      var = c(0),
      residual = c(0),
      zscore = c(0),
      neigh = c(0),
      neighobs = c(0),
      weight = c(0),
      sumWeight = c(0),
      c0arr = c(0)
    )
    cvInfo <- list()
  }
  #

  if (wret) {
    weight <- matrix(0, nrow = npred, ncol = nobs)
  }
  if (interactive() && debug.level == 1 && length(sel) > 1) {
    pb <- txtProgressBar(1, length(sel), style = 3)
  }

  if (debug.level >= 1) {
    print(paste(
      ifelse(cv, "cross-validating", "interpolating "),
      length(sel),
      "areas"
    ))
  }

  for (inn in seq_along(sel)) {
    inew <- sel[inn]
    if (debug.level > 1) {
      print("\n")
    }
    #  for (inew in 1:20) {
    if (interactive() && debug.level == 1 && length(sel) > 1) {
      setTxtProgressBar(pb, inn)
    }

    if (cv) {
      if (debug.level > 1) {
        print(paste(
          "Cross-validating location",
          inew,
          " out of ",
          npred,
          " observation locations"
        ))
      }
      if (debug.level > 1) print(data.frame(observations)[inew, ])
      #      if (cv == inew && inew > 1) browser()
    } else {
      if (debug.level > 1) {
        print(paste(
          "Predicting location ",
          inew,
          " out of ",
          npred,
          " prediction locations"
        ))
      }
      if (debug.level > 1) {
        print(data.frame(predictionLocations)[inew, ])
      }
    }
    newcor <- newcors[inew, ]

    ret <- rkrige(
      data.frame(observations),
      obs0,
      obscors,
      newcor,
      varMatObs,
      varMatPredObs[, inew],
      nmax,
      inew,
      cv,
      unc0,
      mdist,
      maxdist,
      singMat,
      varInv,
      singularSolve,
      wlim,
      debug.level,
      wlimMethod,
      BLUE = BLUE,
      trendObs = trendObs,
      trendPred = trendPred[inew, ]
    )

    predictions$var1.pred[inew] <- ret$pred[1]
    predictions$var1.var[inew] <- ret$pred[2]
    predictions$sumWeights[inew] <- ret$pred[3]

    nneigh <- ret$nneigh
    lambda <- ret$lambda
    neigh <- ret$neigh
    if (wret && !is.na(nneigh)) {
      weight[inew, neigh] <- lambda[1:nneigh]
    }
    obs <- ret$obs
    unc <- ret$unc
    c0arr <- ret$c0arr
    if (debug.level > 1) {
      distm <- utop_dists_n1(obscors, newcor)[neigh]
      lobs <- data.frame(observations)[neigh, ]
      lobs <- rbind(lobs, mu = rep(0, (dim(lobs)[2])))
      lobs <- cbind(
        lobs,
        data.frame(
          id = c(neigh, 0),
          edist = c(distm, 0),
          lambda = lambda[1:(nneigh + 1)],
          c0 = c0arr[1:(nneigh + 1)],
          obs = c(obs, 1),
          unc = c(unc, 0),
          lambda_times_obs = lambda[1:(nneigh + 1)] * c(obs, 0)
        )
      )
      print("neighbours")
      print(lobs, 3)
      print("covariance matrix ")
      print(varMat, 3)
    }
    if (cv) {
      predictions$residual[inew] <- observations[[depVar]][inew] -
        predictions$var1.pred[inew]
      predictions$zscore[inew] <- predictions$residual[inew] /
        sqrt(predictions$var1.var[inew])
      cvNew <- data.frame(
        predictions[inew, ],
        c0arr = lambda[nneigh + 1],
        neigh = 0
      )
      for (jnew in 1:nneigh) {
        cvNew <- rbind(
          cvNew,
          data.frame(
            var1.pred = 0,
            var1.var = unc[neigh[jnew]],
            sumWeights = lambda[jnew],
            observed = data.frame(observations)[[depVar]][neigh[jnew]],
            residual = 0,
            zscore = 0,
            c0arr = c0arr[jnew],
            neigh = neigh[jnew]
          )
        )
      }
      if (inew == 1) {
        cvInfo <- cvNew
      } else {
        cvInfo <- rbind(cvInfo, cvNew)
      }
      if (debug.level > 1) {
        print("prediction")
        print(c(data.frame(predictionLocations)[inew, ], predictions[inew, ]))
      }
    }
    if (singularSolve && !is.null(ret$removed)) {
      if (!is.null(removed)) {
        removed <- rbind(removed, ret$removed)
      } else {
        removed <- ret$removed
      }
    }
  }
  if (interactive() && debug.level == 1 && length(sel) > 1) {
    close(pb)
  }
  predictions <- cbind(predictionLocations, predictions)
  if (cv) {
    predictions$observed <- observations[[depVar]]
  }
  ret <- list(predictions = predictions)
  if (wret) {
    ret$weight <- weight
  }
  if (singularSolve && !is.null(removed)) {
    ret$removed <- removed
  }
  if (cv) {
    ret$cvInfo <- cvInfo
  }
  ret
}
