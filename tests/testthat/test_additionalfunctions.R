test_that("additionalfunctions", {
  x <- cbind(x1 = 3 , x2 = c(4:1, 2:5) - 5)
  expect_equal(mse(x), 7)
  expect_equal(amse(x), 6.625)
})
