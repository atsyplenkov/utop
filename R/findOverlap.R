#' @noRd
findOverlap <- function(areas1, areas2, debug.level = 1, ...) {
  dots <- list(...)
  olim <- if ("olim" %in% names(dots)) dots$olim else 1e-4
  params <- coerce_utop_params(...)
  partialOverlap <- params@partial_overlap

  areas1 <- utop_add_area(areas1)
  if (missing(areas2)) {
    areas2 <- areas1
    sym <- TRUE
  } else {
    areas2 <- utop_add_area(areas2)
    sym <- FALSE
  }

  ndim <- nrow(areas1)
  mdim <- nrow(areas2)
  overlap <- matrix(0, nrow = ndim, ncol = mdim)
  iareas <- utop_area(areas1)
  jareas <- utop_area(areas2)

  if (partialOverlap) {
    for (ia in seq_len(ndim)) {
      ifi <- if (sym) ia else 1
      for (ib in ifi:mdim) {
        inter <- suppressWarnings(sf::st_intersection(
          sf::st_geometry(areas1[ia, ]),
          sf::st_geometry(areas2[ib, ])
        ))
        aover <- if (length(inter) == 0) {
          0
        } else {
          sum(as.numeric(sf::st_area(inter)))
        }
        if (aover / min(iareas[ia], jareas[ib]) > olim) {
          overlap[ia, ib] <- aover
        }
        if (sym) {
          overlap[ib, ia] <- overlap[ia, ib]
        }
      }
    }
    return(overlap)
  }

  pts1 <- suppressWarnings(sf::st_point_on_surface(sf::st_geometry(areas1)))
  pts2 <- suppressWarnings(sf::st_point_on_surface(sf::st_geometry(areas2)))
  geom1 <- sf::st_geometry(areas1)
  geom2 <- sf::st_geometry(areas2)
  nnover <- 0

  for (ia in seq_len(ndim - sym)) {
    t1 <- proc.time()[[3]]
    a1 <- iareas[ia]
    pt1 <- pts1[ia]
    ifi <- ifelse(sym, ia + 1, 1)
    t2 <- proc.time()[[3]]
    for (ib in ifi:mdim) {
      a2 <- jareas[ib]
      if (max(unlist(commonArea(geom1[ia], geom2[ib]))) > 0.1) {
        if (a2 < a1) {
          nover <- sf::st_intersects(pts2[ib], geom1[ia], sparse = FALSE)[1]
        } else {
          nover <- sf::st_intersects(pt1, geom2[ib], sparse = FALSE)[1]
        }
        if (isTRUE(nover)) {
          overlap[ia, ib] <- min(a1, a2)
        }
        nnover <- nnover + 1
        if (sym) {
          overlap[ib, ia] <- overlap[ia, ib]
        }
      }
    }
    t3 <- proc.time()[[3]]
    if (debug.level > 1) {
      print(paste(ia, round(sqrt(a1), 2), round(t2 - t1, 3), round(t3 - t2, 3)))
    }
  }
  if (debug.level > 1) {
    print(paste("nnover", nnover))
  }
  if (sym) {
    diag(overlap) <- iareas
  }
  overlap
}


#' @noRd
find_overlap <- function(
  areas1,
  areas2,
  partial_overlap = TRUE,
  debug.level = 1,
  ...
) {
  findOverlap(
    areas1,
    areas2,
    debug.level = debug.level,
    partialOverlap = partial_overlap,
    ...
  )
}


#' @noRd
findVarioOverlap <- function(vario) {
  overlap <- function(a1, a2, dist) {
    ad <- sqrt(a1) / 2
    ad[2] <- sqrt(a2) / 2
    if (ad[1] + ad[2] > dist) {
      bbox1 <- c(xmin = -ad[1], ymin = -ad[1], xmax = ad[1], ymax = ad[1])
      bbox2 <- c(
        xmin = dist - ad[2],
        ymin = -ad[2],
        xmax = dist + ad[2],
        ymax = ad[2]
      )
      cArea <- commonArea(bbox1, bbox2)
    } else {
      cArea <- list(0, 0)
    }
    cArea[[1]] * a1
  }
  mapply(FUN = overlap, vario$a1, vario$a2, vario$dist)
}


#' @noRd
bbArea <- function(bb) {
  xd <- bb[[3]] - bb[[1]]
  yd <- bb[[4]] - bb[[2]]
  abs(xd) * abs(yd)
}

#' @noRd
commonArea <- function(objecti, objectj) {
  bi <- utop_bbox(objecti)
  bj <- utop_bbox(objectj)
  iarea <- bbArea(bi)
  jarea <- bbArea(bj)
  bl <- list()
  for (i in 1:2) {
    bl[[i]] <- max(bi[[i]], bj[[i]])
  }
  for (i in 3:4) {
    bl[[i]] <- min(bi[[i]], bj[[i]])
  }
  if (bl[[3]] >= bl[[1]] && bl[[4]] >= bl[[2]]) {
    larea <- bbArea(bl)
  } else {
    larea <- 0
  }
  list(larea / iarea, larea / jarea)
}
