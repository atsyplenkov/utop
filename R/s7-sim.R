#' @noRd
utop_sim_runtime_params <- function(params, ...) {
  utop_split_runtime(
    params,
    runtime_names = c(
      "nsim",
      "logdist",
      "beta",
      "dump",
      "large_first",
      "replace",
      "debug_level"
    ),
    ...
  )
}

#' @noRd
compute_sim_utop <- function(
  object,
  var_mat_update = FALSE,
  params = NULL,
  ...
) {
  split <- utop_sim_runtime_params(params, ...)
  object@params <- utop_replace_params(
    current = object@params,
    params = split$params,
    observations = object@observations,
    formula = object@formula
  )
  params_obj <- object@params
  runtime <- split$runtime
  nsim <- if ("nsim" %in% names(runtime)) runtime$nsim else 1L
  beta <- if ("beta" %in% names(runtime)) runtime$beta else NA
  large_first <- if ("large_first" %in% names(runtime)) {
    runtime$large_first
  } else {
    TRUE
  }
  replace <- isTRUE(runtime$replace)
  dump <- runtime$dump
  debug_level <- if ("debug_level" %in% names(runtime)) {
    runtime$debug_level
  } else {
    params_obj@debug_level
  }
  nmax <- params_obj@n_max
  cv <- params_obj@cv
  maxdist <- params_obj@max_dist
  wlim <- params_obj@wlim
  wlim_method <- params_obj@wlim_method
  var_clean <- params_obj@var_clean
  variogram_model <- coerce_variogram_model(object@variogram_model)
  if (is.null(variogram_model)) {
    stop("Cannot do simulations without a variogram model", call. = FALSE)
  }

  if (
    length(object@observations) > 0L &&
      (is.null(object@var_mat_obs) || var_mat_update)
  ) {
    object <- do.call(
      utop_var_mat,
      c(
        list(
          object = object,
          var_mat_update = var_mat_update,
          params = object@params
        ),
        split$dots
      )
    )
  }

  if (
    is.null(object@var_mat_pred) ||
      diff(dim(object@var_mat_pred)) != 0 ||
      var_mat_update
  ) {
    if (
      is.null(object@d_pred) &&
        !(params_obj@g_dist_pred && !is.null(object@g_dist_pred))
    ) {
      object <- do.call(
        utop_disc,
        c(list(object = object, params = object@params), split$dots)
      )
    }
    if (params_obj@g_dist_pred && is.null(object@g_dist_pred)) {
      object <- do.call(
        utop_g_dist,
        c(list(object = object, params = object@params), split$dots)
      )
      object@var_mat_pred <- compute_var_mat_matrix(
        object@g_dist_pred,
        variogramModel = variogram_model
      )
    } else {
      object@var_mat_pred <- compute_var_mat_list(
        object@d_pred,
        variogramModel = variogram_model,
        debug.level = debug_level,
        params = object@params
      )
    }
  }

  var_mat_pred_obs <- object@var_mat_pred_obs
  var_mat_obs <- object@var_mat_obs
  var_mat_pred <- object@var_mat_pred
  predictions <- utop_add_area(object@prediction_locations)
  observations <- utop_add_area(object@observations)
  nobs <- if (is.null(observations)) 0L else nrow(observations)
  prediction_locations <- predictions

  if (replace && !("replaceNumber" %in% names(prediction_locations))) {
    stop(
      "Cannot replace observations if prediction_locations ",
      "does not have column with replaceNumber",
      call. = FALSE
    )
  } else if (
    replace &&
      !(nrow(observations) >=
        max(prediction_locations$replaceNumber, na.rm = TRUE) &&
        min(prediction_locations$replaceNumber, na.rm = TRUE) >= 1)
  ) {
    stop(
      "prediction_locations$replaceNumber does not correspond ",
      "with the number of observations",
      call. = FALSE
    )
  }

  sing_mat <- FALSE
  var_inv <- NULL
  for (isim in seq_len(nsim)) {
    predictions$sim <- NA
    if (length(dim(observations)) > 0 && nrow(observations) > 0) {
      obsall <- data.frame(observations)
      obs <- obsall[, as.character(object@formula[[2]])]
      obscors <- utop_centroid_coordinates(observations)
      nobs0 <- nrow(observations)
    } else {
      obsall <- NULL
      obs <- NULL
      obscors <- NULL
      nobs0 <- 0L
    }
    v_pred <- var_mat_pred
    v_obs <- var_mat_obs
    corlines <- NULL
    if (var_clean) {
      vm <- v_obs
      diag(vm) <- 1
      vm[upper.tri(vm)] <- 1
      mins <- apply(vm, MARGIN = 1, FUN = function(x) min(x))
      corlines <- which(mins < 1e-9)
    }

    v_pred_obs <- var_mat_pred_obs
    d_pred <- nrow(predictions)
    ips <- sample(d_pred, d_pred)
    if (large_first) {
      il <- which(ips == order(predictions$area, decreasing = TRUE)[1])
      tmp <- ips[il]
      ips[il] <- ips[1]
      ips[1] <- tmp
    }
    vpo <- seq_len(d_pred)
    if (interactive() && debug_level) {
      pb <- txtProgressBar(1, length(ips), style = 3)
    }
    print(paste0(isim, ". simulation of ", length(ips), " areas"))
    oind <- NULL
    if (params_obj@unc && "unc" %in% names(observations)) {
      unc0 <- observations$unc
    } else {
      unc0 <- array(0, nobs)
    }

    for (ip in seq_along(ips)) {
      inew <- ips[ip]
      in2 <- which(vpo == inew)
      nobs <- length(obs)
      if (interactive() && debug_level) {
        setTxtProgressBar(pb, ip)
      }
      newcor <- utop_centroid_coordinates(prediction_locations[inew, ])
      if (nobs == 0) {
        if (is.na(beta)) {
          stop(
            "No observations found, beta (expected mean) has to be given",
            call. = FALSE
          )
        }
        c0 <- varioEx(
          sqrt(bbArea(utop_bbox(prediction_locations[in2, ]))),
          variogram_model
        )
        inewvar <- var_mat_pred[inew, inew]
        obs <- rnorm(1, beta, c0 - inewvar)
        v_obs <- matrix(inewvar, nrow = 1, ncol = 1)
        v_pred_obs <- v_pred[inew, -inew, drop = FALSE]
        unc0 <- 0
      } else {
        mdist <- sqrt(diff(range(obscors[, 1]))^2 + diff(range(obscors[, 2]))^2)
        wlim0 <- wlim
        repeat {
          wlim0 <- wlim0 / 1.05
          ret <- try(
            rkrige(
              obsall,
              obs,
              obscors,
              newcor,
              v_obs,
              v_pred_obs[, in2, drop = FALSE],
              nmax,
              inew,
              cv,
              unc0,
              mdist,
              maxdist,
              sing_mat,
              var_inv,
              singularSolve = FALSE,
              wlim0,
              debug_level,
              wlim_method,
              simul = TRUE,
              varClean = FALSE,
              corlines = corlines,
              remNeigh = TRUE
            ),
            silent = TRUE
          )
          if (inherits(ret, "try-error")) {
            print(paste("error in simulation of area number", ip))
          }
          if (
            wlim0 < 1.05 || (!inherits(ret, "try-error") && ret$pred[2] > 0)
          ) {
            break
          }
        }
        if (!inherits(ret, "try-error")) {
          pred <- ret$pred
          newval <- rnorm(1, pred[1], sqrt(pred[2]))
          if (
            replace &&
              !is.na(prediction_locations$replaceNumber)[inew] &&
              prediction_locations$replaceNumber[inew] != 0
          ) {
            obs[prediction_locations$replaceNumber[inew]] <- newval
            oind <- c(oind, prediction_locations$replaceNumber[inew])
          } else {
            obs <- c(obs, newval)
            oind <- c(oind, length(obs))
          }
        } else {
          if (
            replace &&
              !is.na(prediction_locations$replaceNumber)[inew] &&
              prediction_locations$replaceNumber[inew] != 0
          ) {
            obs[prediction_locations$replaceNumber[inew]] <- NA
            oind <- c(oind, prediction_locations$replaceNumber[inew])
          } else {
            obs <- c(obs, NA)
            oind <- c(oind, length(obs))
          }
        }
        if (
          replace &&
            !is.na(prediction_locations$replaceNumber)[inew] &&
            prediction_locations$replaceNumber[inew] != 0
        ) {
          unc0[prediction_locations$replaceNumber[inew]] <- 0
          v_pred_obs <- v_pred_obs[, -in2, drop = FALSE]
        } else {
          unc0 <- c(unc0, 0)
          v_obs <- rbind(v_obs, v_pred_obs[, in2])
          v_obs <- cbind(v_obs, c(v_pred_obs[, in2], 0))
          nd <- nrow(v_obs)
          if (var_clean && any(v_obs[nd, seq_len(nd - 1L)] < 1e-9)) {
            corlines <- c(corlines, nd)
          }
          v_pred_obs <- v_pred_obs[, -in2, drop = FALSE]
          v_pred_obs <- rbind(v_pred_obs, v_pred[in2, -in2])
        }
      }
      obscors <- rbind(obscors, newcor)
      v_pred <- v_pred[-in2, -in2, drop = FALSE]
      vpo <- vpo[-in2]
    }
    if (interactive() && debug_level) {
      close(pb)
    }
    if (replace) {
      predictions$sim[ips] <- obs[oind]
    } else {
      predictions$sim[ips] <- obs[(nobs0 + 1L):length(obs)]
    }
    names(predictions)[ncol(predictions)] <- paste0("sim", isim)
    if (!is.null(dump)) {
      save(object, file = dump)
    }
  }
  object@simulations <- predictions
  object
}

#' Simulate with a utop object
#'
#' @param object A [Utop] object.
#' @param ... Method-specific arguments. For [Utop] objects, `var_mat_update`
#'   and `params` control matrix recomputation and parameter updates.
#'
#' @return A [Utop] object with simulations attached.
#' @export
utop_sim <- S7::new_generic(
  name = "utop_sim",
  dispatch_args = "object",
  fun = function(object, ...) {
    S7::S7_dispatch()
  }
)

S7::method(utop_sim, Utop) <- function(
  object,
  var_mat_update = FALSE,
  params = NULL,
  ...
) {
  compute_sim_utop(
    object,
    var_mat_update = var_mat_update,
    params = params,
    ...
  )
}
