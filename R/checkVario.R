#a) plot 3-d varios
#b) plot upscaled point variograms

#' @noRd
errorBar <- function(x, y, upper, lower = upper, length = 0.1, ...) {
  if (
    length(x) != length(y) ||
      length(y) != length(lower) ||
      length(lower) != length(upper)
  ) {
    stop("vectors must be same length")
  }
  arrows(x, y + upper, x, y - lower, angle = 90, code = 3, length = length, ...)
}

#' @noRd
compute_check_variogram <- function(
  object,
  acor = 1,
  log = "xy",
  cloud = FALSE,
  gDist = TRUE,
  acomp = NULL,
  curveSmooth = FALSE,
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
  dots <- list(...)

  askpar <- par("ask")
  if (dev.interactive()) {
    par("ask" = TRUE)
  } else {
    par("ask" = FALSE)
  }
  variogram_model <- object@variogram_model
  variogramModel <- if (
    S7::S7_inherits(variogram_model, UtopVariogramModel)
  ) {
    coerce_variogram_model(variogram_model)
  } else {
    variogram_model
  }
  # variogram preferred as sampleVariogram, if not variogramCloud (if exisiting, NULL otherwise)
  sampleVariogram <- object@variogram
  if (is.null(sampleVariogram)) {
    sampleVariogram <- object@variogram_cloud
  }
  observations <- object@observations
  formulaString <- object@formula
  amul <- params_obj@amul
  varFit <- object@var_fit
  abins <- adfunc(NULL, observations, amul)
  observations$acl <- findInterval(observations$area, abins)
  observations$n <- 1
  obsvar <- aggregate(
    data.frame(observations)[, as.character(formulaString[[2]])],
    by = list(acl = observations$acl),
    FUN = var
  )
  obsvar$area <- aggregate(
    observations$area,
    by = list(acl = observations$acl),
    FUN = mean
  )[, 2] *
    acor
  obsvar$n <- aggregate(
    observations$n,
    by = list(acl = observations$acl),
    FUN = sum
  )[, 2]
  obsvar$n <- obsvar$n / max(obsvar$n) * 20
  obsvar <- obsvar[!is.na(obsvar$x), ]
  plot(
    obsvar$area,
    obsvar$x,
    xlab = "area",
    ylab = "variance",
    cex = sqrt(obsvar$n),
    pch = 16,
    log = log,
    ylim = c(min(obsvar$x) / 5, max(obsvar$x) * 5)
  )

  print(paste("cloud is ", cloud))
  if (cloud || variogram_is_cloud(sampleVariogram)) {
    print("Creating cloud variogram; this might take some time")
    if (!variogram_is_cloud(sampleVariogram)) {
      if (is.null(object@variogram_cloud)) {
        vario_cloud <- utop_variogram(
          observations,
          formula = formulaString,
          params = params_obj,
          cloud = TRUE
        )
        object@variogram_cloud <- utop_tag_variogram_class(
          utop_variogram_data(vario_cloud),
          cloud = TRUE
        )
      }
      clvar <- object@variogram_cloud
    } else {
      clvar <- sampleVariogram
    }
    if (gDist) {
      if (is.null(object@g_dist_obs)) {
        if (is.null(object@d_obs)) {
          object@d_obs <- utop_disc(observations, params = params_obj)
        }
        dObs <- object@d_obs
        object@g_dist_obs <- compute_g_dist_list(dObs, dObs, params = params_obj)
      }
      gdists <- object@g_dist_obs
      gDiag <- diag(gdists)
      clvar$gdist <- 0
      for (ic in 1:dim(clvar)[1]) {
        ia <- clvar$acl1[ic]
        ib <- clvar$acl2[ic]
        clvar$dist[ic] <- gdists[ia, ib] - 0.5 * (gDiag[ia] + gDiag[ib])
      }
    }
    clvar$np <- clvar$ord
    if (!"identify" %in% names(dots) || !dev.interactive()) {
      dots$identify <- FALSE
    }
    cdots <- which(
      names(dots) %in% names(formals(compute_check_variogram_model))
    )
    if (length(cdots) > 0) {
      dots <- dots[-cdots]
    }
    print(plot(clvar, xlab = "distance", unlist(dots)))
  }
  if (!is.null(varFit) && !variogram_is_cloud(sampleVariogram)) {
    gammar <- varFit[, c("np", "gamma", "gammar")]
    gammar$nnp <- sqrt(gammar$np) / max(sqrt(gammar$np)) * 20
    gmax <- max(gammar[, c("gamma", "gammar")])
    gmin <- quantile(c(gammar$gammar, gammar$gamma), 0.05)
    plot(
      gammar ~ gamma,
      gammar,
      xlim = c(ifelse(length(grep("x", log)) > 0, gmin, 0), gmax),
      ylim = c(ifelse(length(grep("x", log)) > 0, gmin, 0), gmax),
      cex = 0,
      xlab = "gamma",
      ylab = "gamma regularized",
      log = log
    )
    abline(0, 1)
  } else if (!is.null(varFit) && variogram_is_cloud(sampleVariogram)) {
    gammar <- varFit[, c("np", "gamma", "gammar")]
    gammar <- gammar[order(gammar$gamma), ]
    ng <- dim(gammar)[1]
    groups <- ifelse(ng > 200, 20, ng / 10)
    npg <- ng / groups
    gammar$group <- c(1:ng) %/% npg
    ngammar <- aggregate(
      list(gamma = gammar$gamma, gammar = gammar$gammar),
      by = list(gammar$group),
      FUN = mean
    )
    ngammar <- cbind(
      ngammar,
      aggregate(
        list(gammav = gammar$gamma, gammarv = gammar$gammar),
        by = list(gammar$group),
        FUN = var
      )
    )
    xmax <- max(c(ngammar$gamma, ngammar$gammar))
    xmin <- quantile(c(ngammar$gammar, ngammar$gamma), 0.05)

    plot(
      gammar ~ gamma,
      ngammar,
      xlab = "regularized gamma",
      ylab = "gamma",
      xlim = c(ifelse(length(grep("x", log)) > 0, xmin, 0), xmax),
      ylim = c(ifelse(length(grep("x", log)) > 0, xmin, 0), xmax),
      cex = 0,
      pch = 16,
      log = log
    )
    errorBar(ngammar$gammar, ngammar$gamma, upper = sqrt(ngammar$gammav))
    abline(0, 1)
  }
  if (is.null(variogramModel)) {
    if (is.null(sampleVariogram)) {
      sampleVariogram <- utop_variogram(observations, params = params_obj)
    }
  } else {
    if (is.null(sampleVariogram)) {
      object@check_vario <- compute_check_variogram_model(
        object@variogram_model,
        observations = object@observations,
        params = params_obj,
        acor = acor,
        log = log,
        curveSmooth = curveSmooth,
        acomp = acomp,
        ...
      )
    } else {
      object@check_vario <- compute_check_variogram_model(
        object@variogram_model,
        sampleVariogram = sampleVariogram,
        observations = object@observations,
        params = params_obj,
        acor = acor,
        log = log,
        curveSmooth = curveSmooth,
        acomp = acomp,
        ...
      )
    }
  }
  par(ask = askpar)
  invisible(object)
}

#' @noRd
compute_check_variogram_model <- function(
  object,
  sampleVariogram = NULL,
  observations = NULL,
  areas = NULL,
  dists = NULL,
  acomp = NULL,
  params = NULL,
  compVars = list(),
  acor = 1,
  log = "xy",
  legx = NULL,
  legy = NULL,
  plotNugg = TRUE,
  curveSmooth = FALSE,
  ...
) {
  if (S7::S7_inherits(object, UtopVariogramModel)) {
    variogramModel <- coerce_variogram_model(object)
    variogram_model_s7 <- object
  } else {
    variogramModel <- object
    variogram_model_s7 <- utop_variogram_model_from_list(object)
  }
  params_obj <- utop_require_params(params, observations = observations)
  askpar <- par("ask")
  if (dev.interactive()) {
    par("ask" = TRUE)
  } else {
    par("ask" = FALSE)
  }

  if (is.null(areas)) {
    areas <- params_obj@amul
  }
  if (is.null(dists)) {
    dists <- params_obj@dmul
  }

  if (length(areas) == 1) {
    areas <- adfunc(sampleVariogram, observations, areas)
  }
  if (length(dists) == 1) {
    dists <- dfunc(sampleVariogram, observations, dists)
  }

  geoms <- list()
  icomb <- 0
  polylist <- list()
  aavg <- areas[1:(length(areas) - 1)]
  dists <- c(0, dists)
  adists <- dists[seq_along(dists)]
  for (iarea in 1:(length(areas) - 1)) {
    area <- mean(c(areas[iarea], areas[iarea + 1]))
    aavg[iarea] <- area
    for (idist in 1:(length(dists))) {
      icomb <- icomb + 1
      ddist <- ifelse(idist == 1, 0, mean(c(dists[idist], dists[idist - 1])))
      if (iarea == 1) {
        adists[idist] <- ddist
      }
      cs <- sqrt(area) / 2
      x1 <- ddist - cs
      x2 <- ddist + cs
      y1 <- -cs
      y2 <- cs
      boun <- cbind(x = c(x1, x2, x2, x1, x1), y = c(y1, y1, y2, y2, y1))
      geoms[[icomb]] <- sf::st_polygon(list(boun))
    }
  }

  polys <- sf::st_sf(geometry = sf::st_sfc(geoms))
  vmats <- list()
  iplot <- 0
  na <- length(areas)
  if (is.null(acomp) || length(acomp) == 1) {
    if (is.null(acomp)) {
      acomp <- 5
    }
    if (!is.null(sampleVariogram) && !variogram_is_cloud(sampleVariogram)) {
      samp <- aggregate(
        sampleVariogram$np,
        by = list(acl1 = sampleVariogram$acl1, acl2 = sampleVariogram$acl2),
        sum
      )
      if (acomp > dim(samp)[1]) {
        acomp <- dim(samp)[1]
      }
      acomp <- samp[order(samp$x, decreasing = TRUE)[1:acomp], 1:2]
    } else {
      aacomp <- expand.grid(a1 = c(1:(na - 1)), a2 = c(1:(na - 1)))
      aacomp <- aacomp[aacomp$a1 >= aacomp$a2, ]
      if (acomp > dim(aacomp)[1]) {
        acomp <- dim(aacomp)[1]
      }
      acomp <- aacomp[sample(dim(aacomp)[1], acomp), ]
    }
  } else {
    acomp <- acomp[acomp$acl1 < length(areas) & acomp$acl2 < length(areas), ]
  }
  vmats <- matrix(0, ncol = length(dists), nrow = dim(acomp)[1])
  for (iplot in 1:dim(acomp)[1]) {
    i1 <- acomp[iplot, 2]
    i2 <- acomp[iplot, 1]
    ld <- length(adists)
    poly1 <- polys[
      unique(c((i1 - 1) * ld + 1, ((i2 - 1) * ld + 1):(i2 * ld))),
    ]
    poly1$obs <- seq_len(nrow(poly1))
    nadists <- adists
    if (i1 != i2) {
      nadists <- c(0, nadists)
    }
    lobject <- utop_object(poly1, params = params_obj, formula = obs ~ 1)
    lobject@variogram_model <- variogram_model_s7
    overlapObs <- findVarioOverlap(data.frame(
      a1 = utop_area(poly1[1, ]),
      a2 = utop_area(poly1[2, ]),
      dist = nadists
    ))
    lobject@overlap_obs <- t(matrix(
      rep(overlapObs, ld + (i1 != i2)),
      ncol = ld + (i1 != i2)
    ))
    params_cv <- update_utop_params(params_obj, list(cv = TRUE))
    vmat <- utop_var_mat(lobject, params = params_cv)@var_mat_obs
    #  vmat = varMat(poly1,variogramModel = variogramModel, params = params)
    #  vmat = vmat-diag(vmat)
    #    vmats[iplot,] = vmat[1,2:ld]
    #
    if (i1 == i2) {
      vmats[iplot, 2:ld] <- vmat[1, 2:ld]
    } else {
      vmats[iplot, ] <- vmat[1, 2:(ld + 1)]
    }
  }

  if (variogram_is_cloud(sampleVariogram)) {
    xmin <- min(sampleVariogram$dist) / 1.3
  } else {
    xmin <- min(sampleVariogram$dist[sampleVariogram$np > 2] / 1.3)
  }
  xmax <- max(adists)

  pdists <- 10^seq(log10(xmin), log10(xmax), length.out = 100)
  pvar <- apply(
    as.matrix(pdists),
    1,
    varioEx,
    variogramModel = variogramModel
  ) +
    ifelse(plotNugg, nuggEx(1, variogramModel) * acor, 0)
  ymin <- max(min(vmats[vmats > 0]), min(sampleVariogram$gamma))
  ymax <- max(pvar)

  if (acor != 1) {
    xTicks <- axTicks(
      1,
      c(xmin, xmax, 3),
      usr = c(log10(xmin), log10(xmax)),
      log = TRUE,
      nintLog = Inf
    )
    xlabs <- xTicks * sqrt(acor)
  } else {
    xTicks <- NULL
    xlabs <- TRUE
  }
  plot(
    pdists,
    pvar,
    ylim = c(ymin, ymax),
    xlim = c(xmin, xmax),
    log = log,
    type = "l",
    col = "black",
    lwd = 2,
    ylab = "gamma",
    xlab = "distance",
    xaxt = "n"
  )
  axis(1, at = xTicks, labels = xlabs)

  legende <- list(text = "point", col = c("black"), lty = c(1), pch = 16)

  bcols <- c(
    "red",
    "blue",
    "green",
    "orange",
    "brown",
    "violet",
    "yellow",
    "salmon"
  )
  cols1 <- bcols[seq_along(areas)]
  cols2 <- bcols[1:dim(acomp)[1]]
  for (iplot in 1:dim(acomp)[1]) {
    i1 <- acomp[iplot, 2]
    i2 <- acomp[iplot, 1]
    ld <- length(adists)
    if (i1 == i2) {
      lt <- 1
      lcol <- cols1[i1]
    } else {
      lt <- 2
      lcol <- cols2[iplot]
    }

    xx <- adists
    yy <- vmats[iplot, 1:ld]
    if (curveSmooth) {
      if (is.numeric(curveSmooth)) {
        df <- curveSmooth
      } else {
        df <- length(adists) - 3
      }
      xx <- sort(c(xx, seq(min(xx), max(xx), length.out = 1000)))
      yy <- predict(smooth.spline(adists, yy, df = df), xx)$y
    }
    lines(xx, yy, lty = lt, lwd = 2, col = lcol)
    legende$text <- c(
      legende$text,
      paste(aavg[i1] * acor, "vs", aavg[i2] * acor)
    )
    legende$col <- c(legende$col, lcol)
    legende$lty <- c(legende$lty, lt)
    if (!is.null(sampleVariogram) && !variogram_is_cloud(sampleVariogram)) {
      ppts <- sampleVariogram[
        sampleVariogram$acl2 == i1 & sampleVariogram$acl1 == i2,
      ]
      lpch <- 16 + lt
      np <- 0 # Dummy variable to avoid check warning
      points(
        gamma ~ dist,
        ppts,
        col = lcol,
        pch = lpch,
        cex = sqrt(sqrt(np / max(sampleVariogram$np) * 10))
      )
      legende$pch <- c(legende$pch, lpch)
    }
  }
  if (length(compVars) > 0) {
    for (ic in seq_along(compVars)) {
      cvar <- compVars[ic]
      xx <- adists
      if (curveSmooth) {
        xx <- sort(c(xx, seq(min(xx), max(xx), length.out = 1000)))
      }
      clines <- gstat::variogramLine(cvar[[1]], dist_vector = xx)
      lines(clines, lty = 3, lwd = 2, col = cols2[ic])
      legende$text <- c(legende$text, names(cvar))
      legende$col <- c(legende$col, cols2[ic])
      legende$lty <- c(legende$lty, 3)
      legende$pch <- c(legende$pch, 16)
    }
  }
  if (is.null(legx)) {
    legx <- ifelse(
      length(grep("x", log)) > 0,
      max(adists) / log(xmax / xmin, 5),
      max(adists) * 0.7
    )
  }
  if (is.null(legy)) {
    legy <- ifelse(
      length(grep("y", log)) > 0,
      sqrt(ymin * ymax / 1.5),
      ymax * 0.35
    )
  }
  warn <- options("warn")
  options(warn = -1)
  legend(
    legx,
    legy,
    legende$text,
    col = legende$col,
    lty = legende$lty,
    lwd = rep(2, length(legende$pch)),
    pch = legende$pch,
    merge = TRUE
  )

  checkVarioRes <- list(vmats = rbind(vmats, pvar), acomp = acomp)
  par(ask = askpar)
  options(warn = warn$warn)
  invisible(checkVarioRes)
}
