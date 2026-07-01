test_that("UtopParams validators reject invalid states", {
  expect_error(UtopParams(model = c("Ex1", "Sph")), "character scalar")
  expect_error(UtopParams(model = "Bogus"), "not implemented")
  expect_error(UtopParams(r_resol = -1), "positive number")
  expect_error(UtopParams(cloud = c(TRUE, FALSE)), "logical scalar")
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
