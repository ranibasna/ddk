test_that("testing the additionalfunctions", {
  x2 = c(4:1, 2:5)
  x = cbind(x1 = 3 , x2 = c(4:1, 2:5) - 5)
  expect_equal(mse(x2), 1.5)
  expect_equal(amse(x), 6.625)
  f1 <- rbetafda(n = 10, seed = 11)
  expect_equal(as.numeric(format(signif(amse(f1)), nsmall = 1)), 0.109381)
})
