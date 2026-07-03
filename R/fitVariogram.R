#' @noRd
objfunc <- function(
  pars,
  varioIn,
  gDistEst = FALSE,
  dists,
  aOver = NULL,
  model = "Ex1",
  bu,
  bl,
  fit.method = 8,
  debug.level = 0,
  resol = 5,
  nd = 100,
  last = FALSE,
  cloud = NULL,
  ...
) {
  # Debug = 0 means no output
  # Debug = 1 means output for every nd iteration
  # Debug = 2 means output for every iteration
  # Debug = 3 means output for every element in variogram

  variogramModel <- list(model = model, params = pars)
  errSum <- 0
  vario <- data.matrix(varioIn)
  inp <- which(names(varioIn) == "np")
  iacl1 <- which(names(varioIn) == "acl1")
  iacl2 <- which(names(varioIn) == "acl2")
  idist <- which(names(varioIn) == "dist")
  igamma <- which(names(varioIn) == "gamma")
  ia1 <- which(names(varioIn) == "a1")
  ia2 <- which(names(varioIn) == "a2")
  if (is.null(cloud)) {
    cloud <- variogram_is_cloud(varioIn)
  }
  if (gDistEst) {
    #dists is here the n*n matrix gDistObs
    if (dim(dists)[1] == dim(dists)[2]) {
      gd1 <- mapply(
        FUN = function(i, dists) dists[i, i],
        vario[, iacl1],
        MoreArgs = list(dists)
      )
      gd2 <- mapply(
        FUN = function(i, dists) dists[i, i],
        vario[, iacl2],
        MoreArgs = list(dists)
      )
      gb <- mapply(
        FUN = function(i, j, dists) dists[i, j],
        vario[, iacl1],
        vario[, iacl2],
        MoreArgs = list(dists)
      )
      gamma1 <- mapply(
        FUN = varioEx,
        gd1,
        MoreArgs = list(variogramModel = variogramModel)
      )
      gamma2 <- mapply(
        FUN = varioEx,
        gd2,
        MoreArgs = list(variogramModel = variogramModel)
      )
      gammab <- mapply(
        FUN = varioEx,
        gb,
        MoreArgs = list(variogramModel = variogramModel)
      )
      gammar <- gammab - 0.5 * (gamma1 + gamma2)
      if (!is.null(aOver)) {
        farea <- vario[, ia1]
        sarea <- vario[, ia2]
        carea <- mapply(
          FUN = function(i, j, aOver) aOver[i, j],
          vario[, iacl1],
          vario[, iacl2],
          MoreArgs = list(aOver)
        )
        nugg <- mapply(
          FUN = nuggEx,
          (1 / farea + 1 / sarea - 2 * carea / (farea * sarea)) / 2,
          MoreArgs = list(variogramModel = variogramModel)
        )
        gammar <- gammar + nugg
      }
      #         rnugg = (amp/farea+amp/sarea-2.*amp*aov/(farea*sarea))/2.
      #        amp*(1/farea+1/sarea-2*aov/(farea+sarea))/2
    } else {
      # Binned variograms - using geostatistical distance
      # dists is here the n*3 matrix with distances within and between hypothetical areas
      gamma1 <- mapply(
        FUN = varioEx,
        dists[, 1],
        MoreArgs = list(variogramModel = variogramModel)
      )
      gamma2 <- mapply(
        FUN = varioEx,
        dists[, 2],
        MoreArgs = list(variogramModel = variogramModel)
      )
      gammab <- mapply(
        FUN = varioEx,
        dists[, 3],
        MoreArgs = list(variogramModel = variogramModel)
      )
      gammar <- gammab - 0.5 * (gamma1 + gamma2)
      if (!is.null(aOver)) {
        farea <- vario[, ia1]
        sarea <- vario[, ia2]
        carea <- aOver
        nugg <- mapply(
          FUN = nuggEx,
          (1 / farea + 1 / sarea - 2 * carea / (farea * sarea)) / 2,
          MoreArgs = list(variogramModel = variogramModel)
        )
        gammar <- gammar + nugg
      }
    }
  } else {
    # dists is here dObs
    if (cloud) {
      ar1 <- mapply(
        FUN = function(i, dObs) dObs[[i]],
        vario[, iacl1],
        MoreArgs = list(dObs = dists)
      )
      ar2 <- mapply(
        FUN = function(i, dObs) dObs[[i]],
        vario[, iacl2],
        MoreArgs = list(dObs = dists)
      )
      gammar <- mapply(
        FUN = vred,
        a1 = ar1,
        a2 = ar2,
        MoreArgs = list(vredTyp = "ind", variogramModel = variogramModel)
      )
      if (!is.null(aOver)) {
        farea <- vario[, ia1]
        sarea <- vario[, ia2]
        carea <- mapply(
          FUN = function(i, j, aOver) aOver[i, j],
          vario[, iacl1],
          vario[, iacl2],
          MoreArgs = list(aOver)
        )
        nugg <- mapply(
          FUN = nuggEx,
          (1 / farea + 1 / sarea - 2 * carea / (farea * sarea)) / 2,
          MoreArgs = list(variogramModel = variogramModel)
        )
        gammar <- gammar + nugg
      }
    } else {
      # Binned variogram, not geostatistical distance
      gammar <- mapply(
        FUN = vred,
        a1 = vario[, ia1],
        a2 = vario[, ia2],
        dist = vario[, idist],
        MoreArgs = list(
          vredTyp = "hyp",
          variogramModel = variogramModel,
          resol = resol
        )
      )
      if (!is.null(aOver)) {
        farea <- vario[, ia1]
        sarea <- vario[, ia2]
        carea <- aOver
        nugg <- mapply(
          FUN = nuggEx,
          (1 / farea + 1 / sarea - 2 * carea / (farea * sarea)) / 2,
          MoreArgs = list(variogramModel = variogramModel)
        )
        gammar <- gammar + nugg
      }
    }
  }
  err <- mapply(
    FUN = goFit,
    vario[, igamma],
    gammar,
    dist = vario[, idist],
    np = vario[, inp],
    MoreArgs = list(fit.method)
  )
  errSum <- sum(err) / sum(vario[, inp])
  if (debug.level > 1) {
    print(paste(paste(round(pars, 4), collapse = " "), "errSum = ", errSum))
  }
  if (last) {
    return(list(
      varFit = data.frame(vario, vario[, igamma], gammar = gammar, err = err),
      errSum = errSum
    ))
  } else {
    return(errSum)
  }
}


#' @noRd
goFit <- function(gobs, gest, dist, np, fit.method = 8) {
  if (fit.method == 1) {
    ww <- np
  } else if (fit.method == 2) {
    ww <- np * (gobs / gest - 1)^2
    return(ww)
  } else if (fit.method == 6) {
    ww <- 1
  } else if (fit.method == 7) {
    ww <- 1 / dist^2
  } else if (fit.method == 8) {
    ww <- np * (gobs / gest - 1)^2
    return(ww)
  } else if (fit.method == 9) {
    ww <- np * (min((gobs / gest - 1)^2, (gest / gobs - 1)^2))
    return(ww)
  } else {
    stop(paste("fit.method", fit.method, "not implemented"))
  }
  ww * (gobs - gest)^2
}


#' @noRd
varioEx <- function(skor, variogramModel) {
  if (S7::S7_inherits(variogramModel, UtopVariogramModel)) {
    model <- variogramModel@model
    params <- variogramModel@params
  } else {
    model <- variogramModel$model
    params <- variogramModel$params
  }
  res <- 0.
  imod <- imodel(model)
  vres <- .Fortran("varioex", res, skor, length(params), params, imod)
  return(vres[[1]])
}

#' @noRd
imodel <- function(model) {
  #     The numbers should match the numbers of the Fortran-function vario
  as.integer(switch(
    model,
    Exp = 1,
    Ex1 = 2,
    Gau = 3,
    Gho = 5,
    Sph = 6,
    Sp1 = 7,
    Fra = 8
  ))
}

#' @noRd
nuggEx <- function(ared, variogramModel) {
  if (S7::S7_inherits(variogramModel, UtopVariogramModel)) {
    params <- variogramModel@params
  } else {
    params <- variogramModel$params
  }
  return(params[3] * ared)
}
