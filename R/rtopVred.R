#' @noRd
vred <- function(
  a1,
  a2 = NULL,
  vredTyp = "hyp",
  variogramModel,
  pdf1 = NULL,
  pdf2 = NULL,
  aover = NULL,
  dist = NULL,
  inner = 0,
  resol = 5
) {
  #  print(asdf)
  #  return(0)
  model <- variogramModel$model
  imod <- imodel(model)
  param <- variogramModel$params
  ci <- 0
  if (vredTyp == "ind") {
    if (is(a1, "list")) {
      a2 <- a1[[2]]
      a1 <- a1[[1]]
    } else if (is.null(a2)) {
      a2 <- a1
    }
    if (utop_is_legacy_sp(a1) || utop_is_legacy_sp(a2)) {
      utop_stop_legacy_sp(if (utop_is_legacy_sp(a1)) a1 else a2)
    }
    if (inherits(a1, "sf") || inherits(a1, "sfc")) {
      a1 <- utop_point_coordinates(a1)
    }
    if (inherits(a2, "sf") || inherits(a2, "sfc")) {
      a2 <- utop_point_coordinates(a2)
    }
    ip1 <- dim(a1)[1]
    ip2 <- dim(a2)[1]
    vreda <- .Fortran(
      "vredind",
      ci,
      ip1,
      ip2,
      a1,
      a2,
      length(param),
      param,
      imod
    )
    #  } else if (vredTyp == "pdf") {
    #    vreda = .Fortran("vredpdf",ci,c1,c2,ip1,ip2,ipb,pdf1,pdf2,pdfb,length(param),param,model)
  } else if (vredTyp == "hyp") {
    vreda <- .Fortran(
      "vredhyp",
      ci,
      a1,
      a2,
      dist,
      length(param),
      param,
      as.integer(resol),
      imod
    )
  }
  ###### Nugget needs to be implemented
  if (!is.null(aover)) {
    nug <- 0
  } else {
    nug <- 0
  }

  vreda[[1]]
}
