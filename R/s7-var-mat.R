#' @noRd
compute_var_mat_utop <- function(
  object,
  var_mat_update = FALSE,
  full_pred = FALSE,
  params = NULL,
  ...
) {
  object@params <- utop_replace_params(
    current = object@params,
    params = params,
    observations = object@observations,
    formula = object@formula
  )
  params_obj <- object@params
  lg_dist_pred <- params_obj@g_dist_pred

  observations <- if (inherits(object@observations, "stars")) {
    utop_stars_support(object@observations)
  } else {
    utop_as_sf(object@observations)
  }
  n_obs <- nrow(observations)

  prediction_locations <- if (!is.null(object@prediction_locations)) {
    if (inherits(object@prediction_locations, "stars")) {
      utop_stars_support(object@prediction_locations)
    } else {
      utop_as_sf(object@prediction_locations)
    }
  } else {
    NULL
  }

  variogram_model <- coerce_variogram_model(object@variogram_model)
  dots <- list(...)
  debug_level <- if ("debug_level" %in% names(dots)) {
    dots$debug_level
  } else {
    params_obj@debug_level
  }

  a_obs <- utop_area(observations)
  obs_comp <- FALSE
  pred_comp <- FALSE

  if (!is.null(object@var_mat_obs) && !var_mat_update) {
    vm <- object@var_mat_obs
    if (
      !identical(attr(vm, "variogramModel", exact = TRUE), variogram_model) ||
        n_obs != dim(vm)[1] ||
        (!is.null(prediction_locations) &&
          !is.null(object@var_mat_pred_obs) &&
          dim(object@var_mat_pred_obs)[2] != nrow(prediction_locations))
    ) {
      var_mat_update <- TRUE
    }
  }

  if (params_obj@cv && !is.null(object@var_mat_obs) && !var_mat_update) {
    return(object)
  }

  if (is.null(object@var_mat_obs) || var_mat_update) {
    if (
      is.null(object@d_obs) && !(lg_dist_pred && !is.null(object@g_dist_obs))
    ) {
      object <- do.call(
        utop_disc,
        c(list(object = object, params = object@params), dots)
      )
    }
    if (lg_dist_pred) {
      g_dist_obs <- if (!is.null(object@g_dist_obs)) {
        object@g_dist_obs
      } else {
        object <- do.call(
          utop_g_dist,
          c(list(object = object, params = object@params), dots)
        )
        object@g_dist_obs
      }
      if (
        !is.null(params_obj@n_clus) &&
          params_obj@n_clus > 1 &&
          n_obs > params_obj@cn_areas &&
          requireNamespace("parallel", quietly = TRUE)
      ) {
        cl <- utop_cluster_impl(
          params_obj@n_clus,
          type = params_obj@clus_type,
          outfile = params_obj@outfile
        )
        var_mat_obs <- matrix(
          unlist(parallel::parLapply(
            cl,
            seq_len(n_obs),
            fun = function(x, gDistObs, variogramModel) {
              mapply(gDistObs[, x], FUN = function(y) {
                varioEx(y, variogramModel)
              })
            },
            gDistObs = g_dist_obs,
            variogramModel = variogram_model
          )),
          nrow = n_obs,
          ncol = n_obs
        )
      } else {
        var_mat_obs <- matrix(
          mapply(
            g_dist_obs,
            FUN = varioEx,
            MoreArgs = list(variogramModel = variogram_model)
          ),
          nrow = n_obs,
          ncol = n_obs
        )
      }
      v_diag_obs <- diag(var_mat_obs)
      for (ia in seq_len(n_obs - 1L)) {
        for (ib in (ia + 1L):n_obs) {
          var_mat_obs[ia, ib] <- var_mat_obs[ia, ib] -
            0.5 * (v_diag_obs[ia] + v_diag_obs[ib])
          var_mat_obs[ib, ia] <- var_mat_obs[ia, ib]
        }
      }
      object@var_mat_obs <- var_mat_obs
    } else {
      object@var_mat_obs <- compute_var_mat_list(
        object@d_obs,
        coor1 = utop_centroid_coordinates(observations),
        variogramModel = variogram_model,
        debug.level = debug_level,
        params = params_obj
      )
    }
    attr(object@var_mat_obs, "variogramModel") <- variogram_model
    obs_comp <- TRUE
  }

  if (
    !params_obj@cv &&
      !is.null(prediction_locations) &&
      (is.null(object@var_mat_pred_obs) || var_mat_update)
  ) {
    var_mat_obs <- object@var_mat_obs
    v_diag_obs <- diag(var_mat_obs)
    n_pred <- nrow(prediction_locations)

    if (
      is.null(object@d_pred) && !(lg_dist_pred && !is.null(object@g_dist_pred))
    ) {
      object <- do.call(
        utop_disc,
        c(list(object = object, params = object@params), dots)
      )
    }

    if (lg_dist_pred) {
      g_dist_pred <- if (!is.null(object@g_dist_pred)) {
        object@g_dist_pred
      } else {
        object <- do.call(
          utop_g_dist,
          c(list(object = object, params = object@params), dots)
        )
        object@g_dist_pred
      }
      g_dist_pred_obs <- if (!is.null(object@g_dist_pred_obs)) {
        object@g_dist_pred_obs
      } else {
        object <- do.call(
          utop_g_dist,
          c(list(object = object, params = object@params), dots)
        )
        object@g_dist_pred_obs
      }

      print("Creating prediction semivariance matrix. This can take some time.")
      var_mat_pred <- matrix(
        mapply(
          FUN = varioEx,
          g_dist_pred,
          MoreArgs = list(variogramModel = variogram_model)
        ),
        nrow = n_pred,
        ncol = 1L
      )

      if (
        !is.null(params_obj@n_clus) &&
          params_obj@n_clus > 1 &&
          n_obs > params_obj@cn_areas &&
          requireNamespace("parallel", quietly = TRUE)
      ) {
        cl <- utop_cluster_impl(
          n_clus = params_obj@n_clus,
          type = params_obj@clus_type,
          outfile = params_obj@outfile
        )
        var_mat_pred_obs <- matrix(
          unlist(parallel::parLapply(
            cl,
            seq_len(n_pred),
            fun = function(x, gDistPredObs, variogramModel) {
              mapply(gDistPredObs[, x], FUN = function(y) {
                varioEx(y, variogramModel)
              })
            },
            gDistPredObs = g_dist_pred_obs,
            variogramModel = variogram_model
          )),
          nrow = n_obs,
          ncol = n_pred
        )
      } else {
        var_mat_pred_obs <- matrix(
          mapply(
            FUN = varioEx,
            g_dist_pred_obs,
            MoreArgs = list(variogramModel = variogram_model)
          ),
          nrow = n_obs,
          ncol = n_pred
        )
      }

      if (
        is.null(dim(var_mat_pred)) ||
          dim(var_mat_pred)[1] != dim(var_mat_pred)[2]
      ) {
        v_diag_pred <- var_mat_pred
      } else {
        v_diag_pred <- diag(var_mat_pred)
      }
      for (ia in seq_len(n_obs)) {
        for (ib in seq_len(n_pred)) {
          var_mat_pred_obs[ia, ib] <- var_mat_pred_obs[ia, ib] -
            0.5 * (v_diag_obs[ia] + v_diag_pred[ib])
        }
      }
      object@var_mat_pred <- var_mat_pred
      object@var_mat_pred_obs <- var_mat_pred_obs
    } else {
      object@var_mat_pred <- compute_var_mat_list(
        object@d_pred,
        coor1 = utop_centroid_coordinates(prediction_locations),
        diag = TRUE,
        variogramModel = variogram_model,
        debug.level = debug_level,
        params = params_obj
      )
      object@var_mat_pred_obs <- compute_var_mat_list(
        object@d_obs,
        object@d_pred,
        coor1 = utop_centroid_coordinates(observations),
        coor2 = utop_centroid_coordinates(prediction_locations),
        variogramModel = variogram_model,
        sub1 = diag(object@var_mat_obs),
        sub2 = object@var_mat_pred,
        debug.level = debug_level,
        params = params_obj
      )
    }
    pred_comp <- TRUE
  }

  if (params_obj@nugget) {
    if (obs_comp) {
      overlap_obs <- if (!is.null(object@overlap_obs)) {
        object@overlap_obs
      } else {
        object@overlap_obs <- find_overlap(
          observations,
          observations,
          partial_overlap = params_obj@partial_overlap
        )
      }
      f_obs <- matrix(rep(a_obs, n_obs), ncol = n_obs)
      s_obs <- t(f_obs)
      nugg_obs <- matrix(
        mapply(
          FUN = nuggEx,
          (1 / f_obs + 1 / s_obs - 2 * overlap_obs / (f_obs * s_obs)) / 2,
          MoreArgs = list(variogramModel = variogram_model)
        ),
        ncol = n_obs
      )
      diag(nugg_obs) <- 0
      object@var_mat_obs <- object@var_mat_obs + nugg_obs
    }
    if (pred_comp) {
      overlap_pred_obs <- if (!is.null(object@overlap_pred_obs)) {
        object@overlap_pred_obs
      } else {
        object@overlap_pred_obs <- find_overlap(
          observations,
          prediction_locations,
          partial_overlap = params_obj@partial_overlap
        )
      }
      a_pred <- utop_area(prediction_locations)
      n_pred <- length(a_pred)
      f_pred_obs <- matrix(rep(a_obs, n_pred), ncol = n_pred)
      s_pred_obs <- t(matrix(rep(a_pred, n_obs), ncol = n_obs))
      nugg_pred_obs <- matrix(
        mapply(
          FUN = nuggEx,
          (1 /
            f_pred_obs +
            1 / s_pred_obs -
            2 * overlap_pred_obs / (f_pred_obs * s_pred_obs)) /
            2,
          MoreArgs = list(variogramModel = variogram_model)
        ),
        ncol = n_pred
      )
      object@var_mat_pred_obs <- object@var_mat_pred_obs + nugg_pred_obs
    }
  }

  if (!is.null(object@var_mat_pred_obs)) {
    attr(object@var_mat_pred_obs, "variogramModel") <- variogram_model
  }
  if (!is.null(object@var_mat_pred)) {
    attr(object@var_mat_pred, "variogramModel") <- variogram_model
  }

  object
}

#' Compute covariance matrices
#'
#' @param object A [Utop] object, spatial object, matrix, or list.
#' @param ... Method-specific arguments. For [Utop] objects, `var_mat_update`,
#'   `full_pred`, and `params` control matrix recomputation.
#'
#' @return Covariance matrices or an updated [Utop] object.
#' @export
utop_var_mat <- S7::new_generic(
  name = "utop_var_mat",
  dispatch_args = "object",
  fun = function(object, ...) {
    S7::S7_dispatch()
  }
)

S7::method(utop_var_mat, Utop) <- function(
  object,
  var_mat_update = FALSE,
  full_pred = FALSE,
  params = NULL,
  ...
) {
  compute_var_mat_utop(
    object,
    var_mat_update = var_mat_update,
    full_pred = full_pred,
    params = params,
    ...
  )
}

S7::method(utop_var_mat, utop_matrix_class) <- function(object, ...) {
  compute_var_mat_matrix(object, ...)
}

S7::method(utop_var_mat, S7::class_list) <- function(object, ...) {
  compute_var_mat_list(object, ...)
}

#' @noRd
utop_var_mat_prepare_args <- function(args) {
  if ("variogram_model" %in% names(args) && is.null(args$variogramModel)) {
    args$variogramModel <- args$variogram_model
    args$variogram_model <- NULL
  }
  if ("overlap_obs" %in% names(args) && is.null(args$overlapObs)) {
    args$overlapObs <- args$overlap_obs
    args$overlap_obs <- NULL
  }
  if ("overlap_pred_obs" %in% names(args) && is.null(args$overlapPredObs)) {
    args$overlapPredObs <- args$overlap_pred_obs
    args$overlap_pred_obs <- NULL
  }

  inline_params <- setdiff(
    intersect(names(args), utop_param_fields()),
    "params"
  )
  if (length(inline_params) > 0L) {
    stop(
      "pass utop parameters via a UtopParams object in `params`, not as individual arguments: ",
      paste(inline_params, collapse = ", "),
      call. = FALSE
    )
  }
  args$params <- utop_require_params(args$params, observations = args$object1)

  if ("debug_level" %in% names(args) && is.null(args$debug.level)) {
    args$debug.level <- args$debug_level
    args$debug_level <- NULL
  }
  args
}

S7::method(utop_var_mat, utop_sf_class) <- function(object, ...) {
  dots <- list(...)
  object2 <- if ("object2" %in% names(dots)) {
    dots$object2
  } else if (length(dots) > 0L && inherits(dots[[1L]], "sf")) {
    dots[[1L]]
  } else {
    NULL
  }
  if (
    length(dots) > 0L &&
      !("object2" %in% names(dots)) &&
      inherits(dots[[1L]], "sf")
  ) {
    dots <- dots[-1L]
  }
  args <- utop_var_mat_prepare_args(c(
    list(object1 = object, object2 = object2),
    dots
  ))
  result <- do.call(compute_var_mat_default, args)
  if (is.matrix(result)) {
    return(result)
  }
  utop_utop_from_state(list(
    observations = object,
    prediction_locations = object2,
    params = args$params,
    var_mat_obs = result$varMatObs,
    var_mat_pred_obs = result$varMatPredObs,
    var_mat_pred = result$varMatPred,
    overlap_obs = result$overlapObs,
    overlap_pred_obs = result$overlapPredObs
  ))
}

S7::method(utop_var_mat, utop_stars_class) <- function(object, ...) {
  dots <- list(...)
  object2 <- if ("object2" %in% names(dots)) {
    dots$object2
  } else if (length(dots) > 0L && inherits(dots[[1L]], "stars")) {
    dots[[1L]]
  } else {
    NULL
  }
  if (
    length(dots) > 0L &&
      !("object2" %in% names(dots)) &&
      inherits(dots[[1L]], "stars")
  ) {
    dots <- dots[-1L]
  }
  args <- utop_var_mat_prepare_args(c(
    list(
      object1 = utop_stars_support(object),
      object2 = if (is.null(object2)) NULL else utop_stars_support(object2)
    ),
    dots
  ))
  result <- do.call(compute_var_mat_default, args)
  if (is.matrix(result)) {
    return(result)
  }
  utop_utop_from_state(list(
    observations = object,
    prediction_locations = object2,
    params = args$params,
    var_mat_obs = result$varMatObs,
    var_mat_pred_obs = result$varMatPredObs,
    var_mat_pred = result$varMatPred,
    overlap_obs = result$overlapObs,
    overlap_pred_obs = result$overlapPredObs
  ))
}
