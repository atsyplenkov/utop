#' @noRd
utop_stars_dim_names <- function(object) {
  if (!inherits(object, "stars")) {
    stop("Expected a stars object", call. = FALSE)
  }
  dn <- names(dim(object))
  if (length(dn) != 2) {
    stop(
      paste(
        "stars spatiotemporal objects must have exactly two dimensions:",
        "space and time"
      ),
      call. = FALSE
    )
  }
  dn
}

#' @noRd
utop_stars_space_dim <- function(object) {
  dims <- stars::st_dimensions(object)
  is_sfc <- vapply(
    dims,
    FUN = function(x) inherits(x$values, "sfc"),
    FUN.VALUE = logical(1)
  )
  if (sum(is_sfc) != 1) {
    stop(
      "stars objects must have exactly one spatial dimension with sfc values",
      call. = FALSE
    )
  }
  names(dims)[is_sfc]
}

#' @noRd
utop_stars_time_dim <- function(object) {
  dims <- stars::st_dimensions(object)
  space_dim <- utop_stars_space_dim(object)
  dim_names <- setdiff(names(dims), space_dim)
  is_time <- vapply(
    dim_names,
    FUN = function(nm) {
      vals <- stars::st_get_dimension_values(object, nm)
      inherits(vals, c("Date", "POSIXt"))
    },
    FUN.VALUE = logical(1)
  )
  if (sum(is_time) != 1) {
    stop(
      "stars objects must have exactly one Date or POSIXct time dimension",
      call. = FALSE
    )
  }
  dim_names[is_time]
}

#' @noRd
utop_check_stars <- function(object) {
  utop_stars_dim_names(object)
  utop_stars_space_dim(object)
  utop_stars_time_dim(object)
  invisible(object)
}

#' @noRd
utop_stars_nspace <- function(object) {
  dim(object)[[utop_stars_space_dim(object)]]
}

#' @noRd
utop_stars_ntime <- function(object) {
  dim(object)[[utop_stars_time_dim(object)]]
}

#' @noRd
utop_stars_time <- function(object) {
  stars::st_get_dimension_values(object, utop_stars_time_dim(object))
}

#' @noRd
utop_stars_geometry <- function(object) {
  dims <- stars::st_dimensions(object)
  dims[[utop_stars_space_dim(object)]]$values
}

#' @noRd
utop_stars_attr_matrix <- function(object, attr) {
  utop_check_stars(object)
  if (!attr %in% names(object)) {
    stop(paste("stars attribute not found:", attr), call. = FALSE)
  }
  arr <- object[[attr]]
  arr_dim <- dim(arr)
  if (is.null(arr_dim)) {
    stop(paste("stars attribute has no dimensions:", attr), call. = FALSE)
  }
  dim_names <- names(arr_dim)
  if (is.null(dim_names)) {
    dim_names <- names(dim(object))
  }
  space_dim <- utop_stars_space_dim(object)
  time_dim <- utop_stars_time_dim(object)
  perm <- match(c(space_dim, time_dim), dim_names)
  if (anyNA(perm)) {
    stop(
      paste("stars attribute does not use the space/time dimensions:", attr),
      call. = FALSE
    )
  }
  arr <- aperm(arr, perm)
  dim(arr) <- c(utop_stars_nspace(object), utop_stars_ntime(object))
  arr
}

#' @noRd
utop_stars_static_vector <- function(object, attr) {
  mat <- utop_stars_attr_matrix(object, attr)
  vals <- apply(mat, 1, function(x) {
    ux <- unique(x[!is.na(x)])
    if (length(ux) == 0) {
      NA
    } else if (length(ux) == 1) {
      ux
    } else {
      NA
    }
  })
  if (all(is.na(vals))) {
    return(NULL)
  }
  vals
}

#' @noRd
utop_stars_support <- function(object) {
  utop_check_stars(object)
  data <- data.frame(row.names = seq_len(utop_stars_nspace(object)))
  for (nm in names(object)) {
    vals <- utop_stars_static_vector(object, nm)
    if (!is.null(vals)) {
      data[[nm]] <- vals
    }
  }
  support <- sf::st_sf(
    data,
    geometry = utop_stars_geometry(object),
    crs = sf::st_crs(utop_stars_geometry(object))
  )
  if (!"area" %in% names(support)) {
    support$area <- as.numeric(sf::st_area(support))
  }
  support
}

#' @noRd
utop_stars_add_area <- function(object) {
  utop_check_stars(object)
  if ("area" %in% names(object)) {
    return(object)
  }
  area <- as.numeric(sf::st_area(utop_stars_geometry(object)))
  mat <- matrix(
    rep(area, times = utop_stars_ntime(object)),
    nrow = utop_stars_nspace(object)
  )
  utop_stars_set_attr_matrix(object, "area", mat)
}

#' @noRd
utop_stars_set_attr_matrix <- function(object, attr, mat) {
  utop_check_stars(object)
  mat <- as.matrix(mat)
  expected_dim <- c(utop_stars_nspace(object), utop_stars_ntime(object))
  if (!identical(dim(mat), expected_dim)) {
    stop(
      paste("Replacement stars attribute has wrong dimensions:", attr),
      call. = FALSE
    )
  }
  dim_names <- names(dim(object))
  space_dim <- utop_stars_space_dim(object)
  time_dim <- utop_stars_time_dim(object)
  arr <- if (identical(dim_names, c(space_dim, time_dim))) {
    mat
  } else if (identical(dim_names, c(time_dim, space_dim))) {
    t(mat)
  } else {
    stop("Unsupported stars dimension order", call. = FALSE)
  }
  dimnames(arr) <- NULL
  names(dim(arr)) <- dim_names
  object[[attr]] <- arr
  object
}

#' @noRd
utop_stars_slice_time <- function(object, time_index) {
  dim_names <- names(dim(object))
  time_dim <- utop_stars_time_dim(object)
  if (identical(dim_names[1], time_dim)) {
    object[, time_index, drop = FALSE]
  } else {
    object[, , time_index, drop = FALSE]
  }
}

#' @noRd
utop_stars_empty_like <- function(object) {
  stars::st_as_stars(list(), dimensions = stars::st_dimensions(object))
}
