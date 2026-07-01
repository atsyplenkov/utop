expect_snapshot_numeric <- function(x, variant, digits = 8) {
  formatted <- format(
    round(as.numeric(x), digits),
    trim = TRUE,
    scientific = FALSE
  )
  testthat::expect_snapshot_value(
    paste(formatted, collapse = "\n"),
    variant = variant,
    cran = FALSE
  )
}
