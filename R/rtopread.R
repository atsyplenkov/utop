# Remarks
# Is it reasonable to have ainfo <<- read.area.info(finfo,...) to make sure that ainfo
# is also available at the top level after being delivered to function read.areas?

#' create data.frame with observations of data with a spatial support
#'
#' readAreaInfo will read a text file with observations and descriptions of
#' data with a spatial support.
#'
#' The function is of particular use when data are not available as
#' shape-files, or when the observations are not part of the shape-files. This
#' function is mainly for compatibility with the former FORTRAN-version. The
#' simplest way to read the data in that case is through
#' \code{\link[sf]{st_read}}. See also \code{\link{utop-package}}.
#'
#' @param fname name of file with areal information
#' @param id name of column with observation id
#' @param iobs name of column with number of observations
#' @param obs name of column with observations
#' @param unc name of column with possible uncertainty of observation
#' @param filenames name of column with filenames of areas if different names
#' than id should be used.
#' @param sep separator in csv-file
#' @param debug.level used for giving additional output
#' @param moreCols name of other column names the user wants included in ainfo
#' @return \code{data.frame} with information about observations and/or
#' predictionLocations.
#' @author Jon Olav Skoien
#' @keywords spatial
#' @export
readAreaInfo <- function(
  fname = "ainfo.txt",
  id = "id",
  iobs = "iobs",
  obs = "obs",
  unc = "unc",
  filenames = "filenames",
  sep = "\t",
  debug.level = 1,
  moreCols = list(NULL)
) {
  # Separate function to read in information about the areas, with possibility to define column names
  # fname = name of file with information
  # inum = area number
  # id = internal identity - the way areas are stored on hard drive
  # iobs = number of observations for each station
  # obs = the actual observation
  # unc = The uncertainty of the observation, standard deviation
  #       This variable is optional
  # filenames = filenames for areas
  # MoreCols = other variables the user wants to pass on to ainfo
  #  cat(paste(fname))
  if (debug.level > 1) {
    print(paste(fname, id, iobs, obs, unc, filenames, debug.level, moreCols))
  }
  ainfot <- read.csv(fname, header = TRUE, sep = sep)
  if (debug.level > 1) {
    print(summary(ainfot))
  }
  ainfo <- data.frame(
    id = ainfot[, names(ainfot) == id],
    iobs = ainfot[, names(ainfot) == iobs],
    obs = ainfot[, names(ainfot) == obs]
  )
  if (unc %in% names(ainfot)) {
    ainfo <- data.frame(ainfo, unc = ainfot[, names(ainfot) == unc])
  }
  if (filenames %in% names(ainfot)) {
    ainfo <- data.frame(ainfo, filenames = ainfot[, names(ainfot) == filenames])
  }
  # Including remaining arguments, if user wants to include
  ncols <- length(moreCols)
  if (!is.null(moreCols[[1]])) {
    for (icol in 1:ncols) {
      col <- moreCols[[icol]]
      if (debug.level > 1) {
        cat(paste("moreCols", icol, col, "\n"))
      }
      ainfo <- data.frame(ainfo, col1 = ainfot[, names(ainfot) == col])
      names(ainfo) <- c(names(ainfo[1:(length(names(ainfo)) - 1)]), col)
      if (debug.level > 1) {
        cat(paste("moreCols2", icol, col, "\n"))
      }
      if (debug.level > 1) cat(paste("names(ainfo)", names(ainfo), "\n"))
    }
  }
  return(ainfo)
}


#' help file for creating \code{sf} objects with observations and/or
#' predictionLocations of data with a spatial support
#'
#' readAreas will read area-files, add observations and convert the result to
#' \code{\link[sf]{sf}}.
#'
#' If \code{object} is a file name, \code{\link{readAreaInfo}} will be called.
#' If it is a data frame with observations and/or predictionLocations, the
#' function will read areal data from files according to the ID associated with
#' from files according to the ID associated with each
#' observation/predictionLocation.
#'
#' The function is of particular use when data are not available as
#' shape-files, or when the observations are not part of the shape-files. This
#' function is mainly for compatibility with the former FORTRAN-version. The
#' simplest way to read the data in that case is through
#' \code{\link[sf]{st_read}}. See also \code{\link{utop-package}}.
#'
#' @param object either name of file with areal information or data frame with
#' observations
#' @param adir directory where the files with areal information are to be found
#' @param ftype type of file, the only type supported currently is "xy",
#' referring to x- and y-coordinates of boundaries
#' @param projection add projection to the object if input is boundary-files
#' @param ... further parameters to be passed to \code{\link{readAreaInfo}}
#' @return The function creates an \code{sf} object of observations and/or
#' predictionLocations, depending on the information given in \code{object}.
#' and/or predictionLocations, depending on the information given in
#' \code{object}.
#' @author Jon Olav Skoien
#' @keywords spatial
#' @export
readAreas <- function(object, adir = ".", ftype = "xy", projection = NA, ...) {
  # ainfo is 1 - ainfo e.g. read by readAreaInfo
  #          2 - name of the file to pass to readAreaInfo. ainfo is in that case delivered as a top level data.frame
  # pdif gives directory to areal information
  # need option to use other separators as well
  # Output of this function is a list consisting of
  #      1 Updated version of ainfo
  #      2 An sf object with polygons defining the borders of the areas
  if (is.character(object)) {
    cat(paste("calling readAreaInfo with filename ", object, "\n"))
    ainfo <- readAreaInfo(object, ...)
  } else {
    ainfo <- object
  }
  cat(paste(names(ainfo)))
  cat(paste("\n"))
  if (sum(names(ainfo) == "filenames") == 1) {
    fnames <- paste(adir, "/", ainfo$filenames, sep = "")
  } else {
    fnames <- ainfo$id
  }
  row.names(ainfo) <- c(1:dim(ainfo)[1])
  if (ftype == "xy") {
    fnames <- paste(adir, "/", fnames, ".xy", sep = "")
    geometries <- vector("list", length(fnames))
    for (i in seq_along(fnames)) {
      cat(paste("reading first polygon", i, length(fnames), "\n"))
      boun <- read.table(fnames[i], header = FALSE)
      names(boun) <- c("x", "y")
      coords <- as.matrix(boun[, c("x", "y")])
      if (!all(coords[1, ] == coords[nrow(coords), ])) {
        coords <- rbind(coords, coords[1, ])
      }
      geometries[[i]] <- sf::st_polygon(list(coords))
      cat(paste(" Finished polygon\n"))
    }
    crs <- if (is.na(projection)) NA else as.character(projection)
    geom <- sf::st_sfc(geometries, crs = crs)
  } else {
    stop(paste("Filetype", ftype, "not recognized"))
  }
  sfobj <- sf::st_sf(ainfo, geometry = geom)
  sfobj$area <- as.numeric(sf::st_area(sfobj))
  cent <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sfobj)))
  sfobj$labx <- cent[, 1]
  sfobj$laby <- cent[, 2]
  sfobj$bdim <- vapply(geometries, function(poly) nrow(poly[[1]]), numeric(1))
  sfobj
}
