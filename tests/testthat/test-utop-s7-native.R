read_file_text <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

r_source_dir <- function() {
  normalizePath(testthat::test_path("..", "..", "R"), mustWork = FALSE)
}

skip_if_no_source <- function() {
  testthat::skip_if_not(
    file.exists(file.path(r_source_dir(), "utop-package.R")),
    "Package source files unavailable"
  )
}

test_that("S7 implementation does not call legacy conversion helpers", {
  skip_if_no_source()
  r_dir <- r_source_dir()
  files <- list.files(r_dir, pattern = "^s7", full.names = TRUE)
  text <- paste(vapply(files, read_file_text, character(1)), collapse = "\n")

  expect_false(grepl(
    paste(
      "utop_to_rtop",
      "utop_from_rtop",
      "utop_legacy_call",
      "utop_legacy_work_list",
      "utop_apply_legacy_work_list",
      "class\\(work\\) <- \"rtop\"",
      sep = "|"
    ),
    text
  ))
})

test_that("legacy reference adapters file is removed", {
  skip_if_no_source()
  expect_false(file.exists(file.path(r_source_dir(), "legacy-reference-s7.R")))
})

test_that("legacy constructor helpers are removed", {
  skip_if_no_source()
  r_dir <- r_source_dir()
  files <- list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE)
  text <- paste(vapply(files, read_file_text, character(1)), collapse = "\n")

  expect_false(grepl("createRtopObject <- function", text))
  expect_false(grepl("getRtopParams <- function", text))
  expect_false(grepl("getRtopDefaultParams <- function", text))
  expect_false(grepl("rtopCluster <- function", text))
})

test_that("plan legacy symbols are absent from R sources", {
  skip_if_no_source()
  r_dir <- r_source_dir()
  files <- list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE)
  text <- paste(vapply(files, read_file_text, character(1)), collapse = "\n")

  forbidden <- c(
    "UseMethod\\(",
    "@exportS3Method",
    "@rawNamespace S3method",
    "utop_to_rtop\\(",
    "utop_from_rtop\\(",
    "utop_legacy_call\\(",
    "utop_params_rtop_format\\(",
    "utop_variogram_model_rtop_format\\(",
    "class\\(object\\) <- \"rtop\""
  )
  for (pattern in forbidden) {
    expect_false(grepl(pattern, text), info = pattern)
  }
})

test_that("ported S7 modules do not call legacy adapters", {
  skip_if_no_source()
  r_dir <- r_source_dir()
  ported <- c(
    "s7-disc.R",
    "s7-g-dist.R",
    "s7-var-mat.R",
    "s7-krige.R",
    "s7-variogram.R",
    "s7-fit-variogram.R",
    "s7-sim.R",
    "s7-check-variogram.R"
  )
  legacy_calls <- c(
    "utop_params_rtop_format\\(",
    "utop_variogram_model_rtop_format\\(",
    "utop_legacy_compute_state\\(",
    "getRtopParams\\(",
    "createRtopObject\\("
  )
  for (file in ported) {
    text <- read_file_text(file.path(r_dir, file))
    for (pattern in legacy_calls) {
      expect_false(grepl(pattern, text), info = paste(file, pattern))
    }
  }
})

test_that("legacy S3 method bodies are removed from compute cores", {
  skip_if_no_source()
  r_dir <- r_source_dir()
  files <- list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE)
  text <- paste(vapply(files, read_file_text, character(1)), collapse = "\n")

  expect_false(grepl("UseMethod\\(", text))
  expect_false(grepl("@rawNamespace S3method", text))
  expect_false(grepl(
    paste(
      "varMat\\.rtop <- function",
      "rtopKrige\\.rtop <- function",
      "rtopVariogram\\.sf <- function",
      "rtopDisc\\.sf <- function",
      "gDist\\.list <- function",
      "rtopSim\\.rtop <- function",
      "rtopFitVariogram\\.rtop <- function",
      "checkVario\\.rtop <- function",
      "updateRtopVariogram\\.rtop <- function",
      "createRtopObject <- function",
      "getRtopParams <- function",
      sep = "|"
    ),
    text
  ))
})
