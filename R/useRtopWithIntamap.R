#' Legacy intamap integration
#'
#' The old S3 intamap registration path has been removed as part of the S7
#' object model migration. A new integration should be designed against the S7
#' API if intamap support is needed again.
#'
#' @noRd
useRtopWithIntamap <- function() {
  stop(
    "intamap integration for legacy rtop S3 objects has been removed",
    call. = FALSE
  )
}
