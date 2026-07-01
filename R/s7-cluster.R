#' @noRd
utop_cluster_impl <- function(
  n_clus,
  ...,
  action = "start",
  type,
  outfile = NULL
) {
  cl <- getOption("utop.cluster")
  if (length(cl) > 0 && (action == "stop" || action == "restart")) {
    parallel::stopCluster(cl)
    options(utop.cluster = NULL)
  }
  if (length(cl) > 0 && action == "start") {
    if (length(list(...)) > 0) {
      parallel::clusterEvalQ(cl, ...)
    }
  } else if (action == "start" || action == "restart") {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("Not able to start cluster, parallel not available", call. = FALSE)
    }
    if (missing(type) || is.null(type)) {
      cl <- parallel::makeCluster(n_clus, outfile = outfile)
    } else {
      cl <- parallel::makeCluster(n_clus, type, outfile = outfile)
    }
    if (length(list(...)) > 0) {
      parallel::clusterEvalQ(cl, ...)
    }
    options(utop.cluster = cl)
  }
  getOption("utop.cluster")
}

#' Manage a utop parallel cluster
#'
#' Parallel clusters are started automatically when `n_clus > 1` in
#' [utop_params()]. This helper exposes manual start, stop, and restart
#' controls.
#'
#' @param n_clus Number of cluster workers.
#' @param action One of `"start"`, `"stop"`, or `"restart"`.
#' @param type Cluster type passed to [parallel::makeCluster()].
#' @param outfile Optional worker log file.
#' @param ... Commands evaluated on each worker via
#'   [parallel::clusterEvalQ()].
#'
#' @return A cluster object or `NULL`.
#' @export
utop_cluster <- function(n_clus, ..., action = "start", type, outfile = NULL) {
  utop_cluster_impl(
    n_clus = n_clus,
    ...,
    action = action,
    type = type,
    outfile = outfile
  )
}
