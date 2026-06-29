#' @noRd
utop_stop_legacy_sp <- function(object, target = "sf") {
  stop(
    paste(
      "Legacy sp/spacetime objects are no longer supported.",
      "Please convert",
      paste(class(object), collapse = "/"),
      "to",
      target,
      "before calling utop."
    ),
    call. = FALSE
  )
}

#' @noRd
utop_is_legacy_sp <- function(object) {
  inherits(object, "Spatial") || inherits(object, "ST")
}

#' @noRd
utop_as_sf <- function(object, role = "object") {
  if (utop_is_legacy_sp(object)) {
    utop_stop_legacy_sp(object)
  }
  if (inherits(object, "stars")) {
    return(utop_stars_support(object))
  }
  if (inherits(object, "sf")) {
    return(object)
  }
  if (inherits(object, "sfc")) {
    return(sf::st_sf(geometry = object))
  }
  stop(
    paste(
      role,
      "must be an sf, sfc, or stars object; got",
      paste(class(object), collapse = "/")
    ),
    call. = FALSE
  )
}

#' @noRd
utop_area <- function(object) {
  object <- utop_as_sf(object)
  if ("area" %in% names(object)) {
    return(as.numeric(object$area))
  }
  if ("AREA" %in% names(object)) {
    return(as.numeric(object$AREA))
  }
  if ("Shape_Area" %in% names(object)) {
    return(as.numeric(object$Shape_Area))
  }
  as.numeric(sf::st_area(object))
}

#' @noRd
utop_add_area <- function(object) {
  if (inherits(object, "stars")) {
    return(utop_stars_add_area(object))
  }
  object <- utop_as_sf(object)
  if (!"area" %in% names(object)) {
    object$area <- as.numeric(sf::st_area(object))
  }
  object
}

#' @noRd
utop_geometry <- function(object) {
  sf::st_geometry(utop_as_sf(object))
}

#' @noRd
utop_bbox <- function(object) {
  if (is.numeric(object) && length(object) == 4) {
    return(object)
  }
  sf::st_bbox(utop_as_sf(object))
}

#' @noRd
utop_centroid_coordinates <- function(object) {
  object <- utop_as_sf(object)
  coords <- suppressWarnings(sf::st_coordinates(sf::st_centroid(
    sf::st_geometry(object)
  )))
  coords[, 1:2, drop = FALSE]
}

#' @noRd
utop_point_coordinates <- function(object) {
  if (utop_is_legacy_sp(object)) {
    utop_stop_legacy_sp(object)
  }
  if (inherits(object, "sf") || inherits(object, "sfc")) {
    coords <- sf::st_coordinates(sf::st_geometry(sf::st_as_sf(object)))
    return(coords[, 1:2, drop = FALSE])
  }
  if (is.matrix(object)) {
    return(object[, 1:2, drop = FALSE])
  }
  if (is.data.frame(object) && all(c("x", "y") %in% names(object))) {
    return(as.matrix(object[, c("x", "y")]))
  }
  stop(
    paste(
      "Cannot extract point coordinates from",
      paste(class(object), collapse = "/")
    ),
    call. = FALSE
  )
}

#' @noRd
utop_dists_n1 <- function(coords, point) {
  coords <- as.matrix(coords[, 1:2, drop = FALSE])
  point <- as.numeric(point[1:2])
  sqrt(rowSums(sweep(coords, 2, point)^2))
}

#' @noRd
utop_dist_matrix <- function(coords1, coords2 = NULL) {
  coords1 <- as.matrix(coords1[, 1:2, drop = FALSE])
  if (is.null(coords2)) {
    return(as.matrix(stats::dist(coords1)))
  }
  coords2 <- as.matrix(coords2[, 1:2, drop = FALSE])
  dx <- outer(coords1[, 1], coords2[, 1], "-")
  dy <- outer(coords1[, 2], coords2[, 2], "-")
  sqrt(dx^2 + dy^2)
}

#' @noRd
utop_support_size <- function(object) {
  if (inherits(object, "stars")) {
    return(utop_stars_nspace(object))
  }
  nrow(utop_as_sf(object))
}

#' @noRd
utop_data_names <- function(object) {
  if (inherits(object, "stars")) {
    names(object)
  } else if (inherits(object, "sf")) {
    setdiff(names(object), attr(object, "sf_column"))
  } else {
    names(object)
  }
}

#' @noRd
utop_default_formula <- function(object, caller = "formulaString") {
  data_names <- utop_data_names(object)
  if ("obs" %in% data_names) {
    "obs ~ 1"
  } else if ("value" %in% data_names) {
    "value ~ 1"
  } else if (length(data_names) == 1) {
    paste(data_names, "~ 1")
  } else {
    stop(
      paste(caller, "is missing and cannot be found from data"),
      call. = FALSE
    )
  }
}

#' @noRd
utop_clone_with_data <- function(object, data) {
  object <- utop_as_sf(object)
  sf::st_sf(data, geometry = sf::st_geometry(object), crs = sf::st_crs(object))
}
