test_that("UtopParams validators reject invalid states", {
  expect_error(UtopParams(model = c("Ex1", "Sph")), "character scalar")
  expect_error(UtopParams(model = "Bogus"), "not implemented")
  expect_error(UtopParams(r_resol = -1), "greater than 0")
  expect_error(UtopParams(cloud = c(TRUE, FALSE)), "logical scalar")
  expect_error(UtopParams(n_max = 0), "greater than 0")
  expect_error(UtopParams(n_clus = 0), "greater than 0")
  expect_error(UtopParams(wlim = 0), "greater than 0")
  expect_error(UtopParams(uk_trend_support = "bogus"), "centroid, block")
  expect_error(UtopParams(rs_type = "bogus"), "rtop, regular, random")
})

test_that("UtopVariogram validators reject invalid data", {
  expect_error(
    UtopVariogram(data = data.frame(gamma = 1)),
    "missing variogram columns"
  )
  expect_error(UtopVariogramModel(params = numeric()), "non-empty numeric")
})

test_that("UtopVariogram accepts empty data frames", {
  empty <- UtopVariogram(data = data.frame())
  expect_true(S7::S7_inherits(empty, UtopVariogram))
})
