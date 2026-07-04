#' @noRd
compute_var_mat_matrix <- function(
  object,
  variogramModel,
  diag = FALSE,
  sub1,
  sub2,
  ...
) {
  ndim <- dim(object)[1]
  mdim <- dim(object)[2]
  varMatrix <- matrix(
    mapply(FUN = varioEx, object, MoreArgs = list(variogramModel)),
    nrow = ndim,
    ncol = mdim
  )
  if (diag) {
    sub1 <- sub2 <- diag(varMatrix)
  }
  if (!missing(sub1) && !missing(sub2)) {
    for (ia in 1:ndim) {
      for (ib in 1:mdim) {
        if (!(diag & ia == ib)) {
          varMatrix[ia, ib] <- varMatrix[ia, ib] - 0.5 * (sub1[ia] + sub2[ib])
        }
      }
    }
  }
  varMatrix
}

#' @noRd
compute_var_mat_default <- function(
  object1,
  object2 = NULL,
  variogramModel,
  overlapObs,
  overlapPredObs,
  params = NULL,
  ...
) {
  params <- utop_require_params(params, observations = object1)
  variogramModel <- coerce_variogram_model(variogramModel)
  d1 <- utop_disc(object1, params = params, ...)
  if (!is.null(object2)) {
    d2 <- utop_disc(object2, params = params, ...)
  }
  if (params@g_dist_pred) {
    gDist1 <- compute_g_dist_list(d1, params = params, ...)
    # calling varMat.matrix
    varMatObs <- compute_var_mat_matrix(
      gDist1,
      variogramModel = variogramModel,
      diag = TRUE,
      ...
    )
  } else {
    varMatObs <- compute_var_mat_list(
      d1,
      variogramModel = variogramModel,
      params = params,
      ...
    )
  }

  if (params@nugget) {
    if (missing(overlapObs)) {
      overlapObs <- find_overlap(
        object1,
        object1,
        partial_overlap = params@partial_overlap,
        debug.level = params@debug_level
      )
    }

    aObs <- utop_area(object1)
    nObs <- length(aObs)
    nuggObs <- nugget_matrix(aObs, aObs, overlapObs, variogramModel)
    diag(nuggObs) <- 0
    varMatObs <- varMatObs + nuggObs
  }

  if (is.null(object2)) {
    attr(varMatObs, "variogramModel") <- variogramModel
    return(varMatObs)
  }

  if (params@g_dist_pred && !is.null(object2)) {
    gDistPred <- compute_g_dist_list(d2, diag = TRUE, params = params, ...)
    # Calling varMat.matrix
    varMatPred <- compute_var_mat_matrix(
      gDistPred,
      variogramModel = variogramModel,
      ...
    )
    gDistPredObs <- compute_g_dist_list(d1, d2, params = params, ...)
    varMatPredObs <- compute_var_mat_matrix(
      gDistPredObs,
      sub1 = diag(varMatObs),
      sub2 = varMatPred,
      variogramModel = variogramModel,
      ...
    )
  } else {
    varMatPred <- compute_var_mat_list(
      d2,
      diag = TRUE,
      params = params,
      variogramModel = variogramModel
    )
    varMatPredObs <- compute_var_mat_list(
      d1,
      d2,
      sub1 = diag(varMatObs),
      sub2 = varMatPred,
      params = params,
      variogramModel = variogramModel,
      ...
    )
  }
  if (params@nugget) {
    if (missing(overlapPredObs)) {
      overlapPredObs <- find_overlap(
        object1,
        object2,
        partial_overlap = params@partial_overlap,
        debug.level = params@debug_level
      )
    }

    aPred <- utop_area(object2)
    nuggPredObs <- nugget_matrix(aObs, aPred, overlapPredObs, variogramModel)
    varMatPredObs <- varMatPredObs + nuggPredObs
  }
  attr(varMatObs, "variogramModel") <- variogramModel
  attr(varMatPredObs, "variogramModel") <- variogramModel
  attr(varMatPred, "variogramModel") <- variogramModel
  return(list(
    varMatObs = varMatObs,
    varMatPred = varMatPred,
    varMatPredObs = varMatPredObs
  ))
}


# object and object2 (as lists) are here discretized areas
# coor1 and coor2 are coordinates of the areas, used for maximum distance
#' @noRd
compute_var_mat_list <- function(
  object,
  object2 = NULL,
  coor1,
  coor2,
  maxdist = Inf,
  variogramModel,
  diag = FALSE,
  sub1,
  sub2,
  params = NULL,
  debug.level = NULL,
  debug_level = NULL,
  ...
) {
  params <- utop_require_params(params)
  if (is.null(debug.level)) {
    debug.level <- if (!is.null(debug_level)) {
      debug_level
    } else {
      params@debug_level
    }
  }
  d1 <- object
  d2 <- object2
  if (is.null(d2)) {
    equal <- TRUE
    d2 <- d1
    if (!missing(coor1)) coor2 <- coor1
  } else {
    equal <- FALSE
  }
  if (diag) {
    lmat <- mapply(
      vred,
      a2 = d1,
      a1 = d1,
      MoreArgs = list(vredTyp = "ind", variogramModel = variogramModel)
    )
    return(lmat)
  }

  ndim <- length(d1)
  mdim <- length(d2)
  varMatrix <- matrix(-999, nrow = ndim, ncol = mdim)

  if (
    !is.null(params@n_clus) &&
      params@n_clus > 1 &&
      length(d1) + length(d2) > params@cn_areas &&
      requireNamespace("parallel")
  ) {
    cl <- utop_cluster_impl(
      params@n_clus,
      type = params@clus_type,
      outfile = params@outfile
    )
    if (missing(coor1)) {
      coor1 <- NULL
    }
    if (missing(coor2)) {
      coor2 <- NULL
    }

    fun <- function(ia, d1, d2, coor1, coor2, equal, maxdist, debug.level) {
      if (debug.level > 10) {
        print(paste("ia", typeof(ia)))
      }
      t1 <- proc.time()[[3]]
      ndim <- length(d1)
      if (equal) {
        mdim <- ndim
      } else {
        mdim <- length(d2)
      }
      a1 <- d1[[ia]]
      a1 <- utop_point_coordinates(a1)

      first <- ifelse(equal, ia, 1)
      lorder <- c(first:mdim)
      if (!is.null(coor1) && !is.null(coor2) && maxdist < Inf) {
        lorder <- lorder[
          utop_dists_n1(coor2[first:mdim, ], coor1[ia, ]) < maxdist
        ]
      }
      if (length(lorder) > 0) {
        if (equal) {
          a2 <- d1[lorder]
        } else {
          a2 <- d2[lorder]
        }

        lmat <- mapply(
          vred,
          a2 = a2,
          MoreArgs = list(
            vredTyp = "ind",
            a1 = a1,
            variogramModel = variogramModel
          )
        )
      } else {
        lmat <- -999
      }
      t2 <- proc.time()[[3]]
      if (debug.level > 0) {
        print(paste(
          "varMat - Finished element ",
          ia,
          " in ",
          round(t2 - t1, 3),
          "PID:",
          Sys.getpid()
        ))
      }
      list(lmat, lorder)
    }
    if (equal) {
      d2 <- NULL
    }
    parallel::clusterExport(
      cl,
      c(
        "d1",
        "d2",
        "coor1",
        "coor2",
        "equal",
        "maxdist",
        "fun",
        "debug.level",
        "variogramModel",
        "utop_point_coordinates",
        "utop_dists_n1",
        "vred"
      ),
      envir = environment()
    )
    vmll <- parallel::clusterApplyLB(
      cl,
      seq_along(d1),
      d1 = d1,
      d2 = d2,
      coor1 = coor1,
      coor2 = coor2,
      equal = equal,
      maxdist = maxdist,
      fun = fun,
      debug.level = debug.level
    )

    for (ia in seq_along(vmll)) {
      lmat <- vmll[[ia]][[1]]
      lorder <- vmll[[ia]][[2]]
      if (!equal) {
        varMatrix[ia, ] <- lmat
      } else {
        varMatrix[ia, lorder] <- lmat
        varMatrix[lorder, ia] <- lmat
      }
    }
  } else {
    t0 <- proc.time()[[3]]
    for (ia in 1:ndim) {
      t1 <- proc.time()[[3]]
      a1 <- d1[[ia]]
      a1 <- utop_point_coordinates(a1)
      first <- ifelse(equal, ia, 1)
      lorder <- c(first:mdim)
      if (!missing(coor1) && !missing(coor2) && maxdist < Inf) {
        lorder <- lorder[
          utop_dists_n1(coor2[first:mdim, ], coor1[ia, ]) < maxdist
        ]
      }
      if (length(lorder) > 0) {
        a2 <- d2[lorder]
        lmat <- mapply(
          vred,
          a2 = a2,
          MoreArgs = list(
            vredTyp = "ind",
            a1 = a1,
            variogramModel = variogramModel
          )
        )
        if (!equal) {
          varMatrix[ia, ] <- lmat
        } else {
          varMatrix[ia, lorder] <- lmat
          varMatrix[lorder, ia] <- lmat
        }
      }
      t2 <- proc.time()[[3]]
      if (debug.level > 0) {
        print(paste(
          "varMat - Finished element ",
          ia,
          " out of ",
          ndim,
          " in ",
          round(t2 - t1, 3),
          "seconds - totally",
          round(t2 - t0),
          " seconds"
        ))
      }
    }
  }
  if (equal && variogramModel$model != "Gho") {
    varMatrix <- regularize_symmetric(varMatrix)
  } else if (!missing(sub1) && !missing(sub2)) {
    varMatrix <- regularize_cross(varMatrix, sub1, sub2)
  }
  attr(varMatrix, "variogramModel") <- variogramModel
  varMatrix
}
