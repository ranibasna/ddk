test_that("testing the additionalfunctions", {
  x2 = c(4:1, 2:5)
  x = cbind(x1 = 3 , x2 = c(4:1, 2:5) - 5)
  expect_equal(mse(x2), 1.5)
  expect_equal(amse(x), 6.625)
  f1 <- rbetafda(n = 10, seed = 11)
  expect_equal(as.numeric(format(signif(amse(f1)), nsmall = 1)), 0.109381)
  # test vector data class
  f_vec <- rbetafda(n = 1, nx = 10, seed = 3)
  expect_equal(amse(f_vec), 0.004756012)
  expect_equal(mse(f_vec), 0.004756012)
})
