test_that("testing the opt_split function", {
  f <- rbetafda(n = 4, nx = 4)
  AMSE = amse(f)
  expect_equal(class(opt_split(f = f, AMSE = AMSE, M = 2)), "list")
  expect_equal(length(opt_split(f = f, AMSE = AMSE, M = 2)), 3)
})
