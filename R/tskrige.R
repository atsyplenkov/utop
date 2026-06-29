#' @export
#' @rdname rtopKrige
rtopKrige.stars <- function(
  object,
  predictionLocations = NULL,
  varMatObs,
  varMatPredObs,
  varMat,
  params = list(),
  formulaString,
  sel,
  olags = NULL,
  plags = NULL,
  lagExact = TRUE,
  ...
) {
  params <- getRtopParams(params, ...)
  cv <- params$cv
  debug.level <- params$debug.level
  if (!cv && !isTRUE(all.equal(is.null(olags), is.null(plags)))) {
    stop(paste(
      "Lag times have to be given for both observations and",
      "predictionLocations, or for none of them"
    ))
  }
  if (!missing(varMat) && missing(varMatObs)) {
    if (is.atomic(varMat)) {
      varMatObs <- varMat
      if (!cv) {
        stop(paste(
          "Not cross-validation, you must provide either varMatObs",
          "and varMatPredObs or varMat with the matrices as elements"
        ))
      }
    } else {
      varMatObs <- varMat$varMatObs
      varMatPredObs <- varMat$varMatPredObs
    }
  }
  if (missing(formulaString)) {
    formulaString <- utop_default_formula(object)
  }
  if (!inherits(formulaString, "formula")) {
    formulaString <- as.formula(formulaString)
  }

  object <- utop_stars_add_area(object)
  if (cv) {
    predictionLocations <- object
    plags <- olags
    varMatPredObs <- varMatObs
  } else if (is.null(predictionLocations)) {
    stop("predictionLocations are required unless cv = TRUE")
  } else {
    predictionLocations <- utop_stars_add_area(predictionLocations)
  }

  obs_support <- utop_stars_support(object)
  pred_support <- if (cv) {
    obs_support
  } else {
    utop_stars_support(predictionLocations)
  }

  tObs <- ukTrendMatrix(formulaString, obs_support, params)
  tPred <- if (cv) {
    tObs
  } else {
    ukTrendMatrix(formulaString, pred_support, params)
  }

  depVar <- as.character(formulaString[[2]])
  obs <- utop_stars_attr_matrix(object, depVar)
  obs2 <- obs
  nspace <- utop_stars_nspace(object)
  ntime <- utop_stars_ntime(object)
  pspace <- utop_stars_nspace(predictionLocations)
  ptime <- utop_stars_ntime(predictionLocations)
  time_match <- match(
    utop_stars_time(predictionLocations),
    utop_stars_time(object)
  )

  pred_mat <- matrix(NA_real_, nrow = pspace, ncol = ptime)
  var_mat <- matrix(NA_real_, nrow = pspace, ncol = ptime)
  yam_mat <- matrix(NA_real_, nrow = pspace, ncol = ptime)

  if (is.null(olags) && !anyNA(obs)) {
    obs_first <- obs_support
    obs_first[[depVar]] <- obs[, 1]
    ret <- rtopKrige.default(
      obs_first,
      pred_support,
      varMatObs,
      varMatPredObs,
      varMat,
      params,
      formulaString,
      wret = TRUE,
      trendObs = tObs,
      trendPred = tPred
    )
    weight <- ret$weight
    wvar <- ret$predictions$var1.var
    for (istat in seq_len(pspace)) {
      lweight <- weight[istat, ]
      for (itime in seq_len(ptime)) {
        otime <- time_match[itime]
        if (is.na(otime)) {
          next
        }
        pred <- sum(lweight * obs[, otime])
        pred_mat[istat, itime] <- pred
        var_mat[istat, itime] <- wvar[istat]
        yam_mat[istat, itime] <- sum(lweight * (obs[, otime] - pred)^2)
      }
    }
  } else if (is.null(olags)) {
    oldind <- NULL
    for (itime in seq_len(ptime)) {
      otime <- time_match[itime]
      if (is.na(otime)) {
        next
      }
      newind <- which(!is.na(obs[, otime]))
      if (length(newind) == 0) {
        next
      }
      if (is.null(oldind) || !isTRUE(all.equal(newind, oldind))) {
        oldind <- newind
        ppq <- obs_support[newind, ]
        ppq[[depVar]] <- obs[newind, otime]
        vmat <- varMatObs[newind, newind]
        vpredobs <- if (cv) {
          NULL
        } else {
          varMatPredObs[newind, , drop = FALSE]
        }
        ret <- rtopKrige.default(
          object = ppq,
          predictionLocations = if (cv) NULL else pred_support,
          varMatObs = vmat,
          varMatPredObs = vpredobs,
          params = params,
          formulaString = formulaString,
          wret = TRUE,
          debug.level = 0,
          trendObs = tObs[newind, , drop = FALSE],
          trendPred = if (cv) NULL else tPred
        )
        weight <- ret$weight
        wvar <- ret$predictions$var1.var
      }
      ob <- obs[newind, otime]
      preds <- weight %*% ob
      if (cv) {
        target <- newind
      } else {
        target <- seq_len(pspace)
      }
      pred_mat[target, itime] <- as.vector(preds)
      var_mat[target, itime] <- wvar
      for (istat in seq_along(target)) {
        diffs <- ob - preds[istat]
        yam_mat[target[istat], itime] <- sum(weight[istat, ] * (diffs^2))
      }
    }
  } else {
    obs_coords <- utop_centroid_coordinates(obs_support)
    pred_coords <- utop_centroid_coordinates(pred_support)
    if (is.null(plags)) {
      stop("plags must be supplied when olags are supplied")
    }
    if (length(olags) != nspace || length(plags) != pspace) {
      stop(paste(
        "olags and plags must match the number of",
        "observation/prediction supports"
      ))
    }
    for (istat in seq_len(pspace)) {
      ppred <- pred_support[istat, ]
      rolags <- olags - plags[istat]
      nbefore <- max(0, -floor(min(rolags)))
      nafter <- ceiling(max(c(rolags, 1)))
      obsb <- cbind(
        matrix(NA_real_, nrow = nspace, ncol = nbefore),
        obs2,
        matrix(NA_real_, nrow = nspace, ncol = nafter)
      )
      obs[,] <- NA_real_
      if (!lagExact) {
        rolags <- round(rolags, 0)
        for (jstat in seq_len(nspace)) {
          obs[jstat, ] <- obsb[
            jstat,
            (nbefore + 1 + rolags[jstat]):(nbefore + ntime + rolags[jstat])
          ]
        }
      } else {
        rdiff <- rolags - floor(rolags)
        for (jstat in seq_len(nspace)) {
          nf <- (nbefore + 1 + floor(rolags[jstat])):(nbefore +
            ntime +
            floor(rolags[jstat]))
          nl <- (nbefore + 1 + ceiling(rolags[jstat])):(nbefore +
            ntime +
            ceiling(rolags[jstat]))
          obs[jstat, ] <- obsb[jstat, nf] *
            (1 - rdiff[jstat]) +
            obsb[jstat, nl] * rdiff[jstat]
        }
      }
      if (cv) {
        vorder <- order(varMatObs[istat, ])
      } else {
        # varMatPredObs stores observations in rows and predictions in columns.
        vorder <- order(varMatPredObs[, istat])
      }
      oldind <- NULL
      for (itime in seq_len(ptime)) {
        jtime <- time_match[itime]
        if (is.na(jtime)) {
          next
        }
        newind <- vorder[vorder %in% which(!is.na(obs[, jtime]))]
        if (cv) {
          newind <- newind[!newind %in% istat]
        }
        if (is.numeric(params$nmax) && length(newind) > params$nmax) {
          newind <- newind[seq_len(params$nmax)]
        }
        if (length(newind) == 0) {
          next
        }
        if (is.null(oldind) || !isTRUE(all.equal(newind, oldind))) {
          oldind <- newind
          ppq <- obs_support[newind, ]
          ppq$intvar <- obs[newind, jtime]
          vmat <- varMatObs[newind, newind]
          vpredobs <- if (cv) {
            varMatObs[newind, istat, drop = FALSE]
          } else {
            varMatPredObs[newind, istat, drop = FALSE]
          }
          ret <- rtopKrige.default(
            object = ppq,
            ppred,
            varMatObs = vmat,
            varMatPredObs = vpredobs,
            params = params,
            formulaString = intvar ~ 1,
            wret = TRUE,
            debug.level = 0,
            cv = FALSE,
            trendObs = tObs[newind, , drop = FALSE],
            trendPred = tPred[istat, , drop = FALSE]
          )
          weight <- ret$weight
          wvar <- ret$predictions$var1.var
        }
        pred <- sum(weight * obs[newind, jtime])
        pred_mat[istat, itime] <- pred
        var_mat[istat, itime] <- wvar
        yam_mat[istat, itime] <- sum(weight * (obs[newind, jtime] - pred)^2)
      }
    }
    if (debug.level > 2) {
      print(list(obs_coords = obs_coords, pred_coords = pred_coords))
    }
  }

  predictions <- predictionLocations
  predictions <- utop_stars_set_attr_matrix(predictions, "var1.pred", pred_mat)
  predictions <- utop_stars_set_attr_matrix(predictions, "var1.var", var_mat)
  predictions <- utop_stars_set_attr_matrix(predictions, "var1.yam", yam_mat)
  list(predictions = predictions)
}
