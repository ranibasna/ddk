test_that("rbetafda", {
  f1 <- rbetafda(n = 10, seed = 11)
  expect_equal(dim(f1), c(10,500))
})
