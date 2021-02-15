test_that("testing the split function", {
  f <- rbetafda(n = 4, nx = 4, seed = 11)
  expect_equal(length(split(f,2)), 2)
  expect_equal(split(f,2), c(0.003040530, 0.002576843))
})
