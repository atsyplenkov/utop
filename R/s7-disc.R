#' @noRd
compute_disc_sf <- function(
  object,
  params,
  bb = sf::st_bbox(object),
  ...
) {
  if (!S7::S7_inherits(params, UtopParams)) {
    params <- coerce_utop_params(params, ...)
  }
  stype <- params@rs_type
  resol <- params@r_resol
  debug.level <- params@debug_level
  if (stype == "random" || stype == "regular") {
    return(lapply(sf::st_geometry(object), FUN = function(pol) {
      sf::st_sample(pol, size = resol, type = stype, offset = c(0.5, 0.5))
    }))
  }
  if (stype != "rtop") {
    stop(paste("Unknown sampling type:", stype), call. = FALSE)
  }

  bbdia <- sqrt(bbArea(bb))
  small <- bbdia / 100
  ires0 <- 1
  nps <- dim(object)[1]
  spp <- vector("list", nps)

  lfun <- function(lpoly, resol, ires0, bbdia, small) {
    if (!is.na(sf::st_crs(lpoly))) {
      lpoly <- sf::st_set_crs(lpoly, NA)
    }
    ba <- sf::st_bbox(lpoly)
    ipts <- resol - 1
    ires <- ires0
    while (ipts < resol) {
      ires <- ires * 2
      xd <- bbdia / (ires)
      if (bbArea(ba) / (xd * xd) > (resol - 2)) {
        x <- seq(bb[[1]] - small, bb[[3]] + small, xd)
        y <- seq(bb[[2]] - small, bb[[4]] + small, xd)
        x <- x[x > ba[[1]] & x < ba[[3]]]
        y <- y[y > ba[[2]] & y < ba[[4]]]
        pts <- expand.grid(x = x, y = y)
        if (dim(pts)[1] >= 1) {
          pts <- sf::st_as_sf(pts, coords = c("x", "y"))
          pts <- pts[sf::st_intersects(lpoly, pts)[[1]], 1]
          ipts <- dim(pts)[1]
        }
      }
    }
    pts
  }

  if (
    !is.null(params@n_clus) &&
      params@n_clus > 1 &&
      dim(object)[1] * params@r_resol / 100 > params@cn_areas
  ) {
    if (!suppressMessages(suppressWarnings(requireNamespace("parallel")))) {
      stop("nclus is > 1, but package parallel is not available", call. = FALSE)
    }
    cl <- utop_cluster_impl(
      params@n_clus,
      type = params@clus_type,
      outfile = params@outfile
    )
    parallel::clusterExport(
      cl,
      c("resol", "ires0", "bbdia", "small"),
      envir = environment()
    )
    spp <- parallel::clusterApply(
      cl,
      sf::st_geometry(object),
      fun = function(x) {
        lfun(x, resol, ires0, bbdia, small)
      }
    )
  } else {
    if (interactive() && debug.level <= 1) {
      pb <- txtProgressBar(1, nps, style = 3)
    }
    print(paste("Sampling points from ", nps, "areas"))
    for (ip in 1:nps) {
      spp[[ip]] <- lfun(
        sf::st_geometry(object)[ip],
        resol,
        ires0,
        bbdia,
        small
      )
      if (debug.level > 1) {
        print(paste(
          "Sampling from area number",
          ip,
          "containing",
          dim(spp[[ip]])[1],
          "points"
        ))
      } else if (interactive()) {
        setTxtProgressBar(pb, ip)
      }
    }
    if (interactive() && debug.level <= 1) {
      close(pb)
    }
    if (debug.level >= 0) {
      print(paste(
        "Sampled on average",
        round(
          mean(unlist(lapply(spp, FUN = function(sppp) dim(sppp)[1]))),
          2
        ),
        "points from",
        nps,
        "areas"
      ))
    }
  }
  spp
}

#' @noRd
compute_disc_stars <- function(object, params, bb = NULL, ...) {
  support <- utop_stars_support(object)
  if (is.null(bb)) {
    bb <- sf::st_bbox(support)
  }
  compute_disc_sf(support, params = params, bb = bb, ...)
}

#' @noRd
compute_disc_variogram <- function(variogram, params, ...) {
  if (!S7::S7_inherits(params, UtopParams)) {
    params <- coerce_utop_params(params, ...)
  }
  data <- utop_variogram_data(variogram)
  resol <- params@h_resol^2
  hstype <- params@hs_type
  mapply(
    disc_bin_areas,
    as.list(data$a1),
    as.list(data$a2),
    as.list(data$dist),
    MoreArgs = list(resol = resol, stype = hstype),
    SIMPLIFY = FALSE
  )
}

#' Discretise spatial supports
#'
#' @param object A [Utop] object, spatial object, or variogram object.
#' @param ... Method-specific arguments. For [Utop] objects, `params` updates
#'   parameters.
#'
#' @return Discretisation points or an updated [Utop] object.
#' @export
utop_disc <- S7::new_generic(
  name = "utop_disc",
  dispatch_args = "object",
  fun = function(object, ...) {
    S7::S7_dispatch()
  }
)

S7::method(utop_disc, Utop) <- function(object, params = list(), ...) {
  object@params <- utop_params(object@params, new_params = params, ...)
  bb <- utop_support_bbox(
    object@observations,
    object@prediction_locations
  )

  object@d_obs <- utop_disc(
    object@observations,
    bb = bb,
    params = object@params,
    ...
  )

  if (inherits(object@observations, "sf")) {
    object@observations$ddim <- vapply(
      object@d_obs,
      FUN = function(are) nrow(utop_point_coordinates(are)),
      FUN.VALUE = numeric(1)
    )
  }

  if (!is.null(object@prediction_locations)) {
    object@d_pred <- utop_disc(
      object@prediction_locations,
      bb = bb,
      params = object@params,
      ...
    )
    if (inherits(object@prediction_locations, "sf")) {
      object@prediction_locations$ddim <- vapply(
        object@d_pred,
        FUN = function(are) nrow(utop_point_coordinates(are)),
        FUN.VALUE = numeric(1)
      )
    }
  }

  object
}

S7::method(utop_disc, UtopVariogram) <- function(object, params = list(), ...) {
  if (!S7::S7_inherits(params, UtopParams)) {
    params <- utop_params(params, ...)
  }
  compute_disc_variogram(object, params = params, ...)
}

S7::method(utop_disc, utop_sf_class) <- function(
  object,
  params = UtopParams(),
  bb = sf::st_bbox(object),
  ...
) {
  if (!S7::S7_inherits(params, UtopParams)) {
    params <- utop_params(params, ...)
  }
  compute_disc_sf(object, params = params, bb = bb, ...)
}

S7::method(utop_disc, utop_stars_class) <- function(
  object,
  params = UtopParams(),
  bb = NULL,
  ...
) {
  if (!S7::S7_inherits(params, UtopParams)) {
    params <- utop_params(params, ...)
  }
  compute_disc_stars(object, params = params, bb = bb, ...)
}
