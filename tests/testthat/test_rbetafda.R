test_that("testing the rbetafda", {
  f1 <- rbetafda(n = 10, seed = 11)
  expect_equal(dim(f1), c(10,500))
  expect_equal(class(f1), c("matrix","array"))
  f2 <- rbetafda(n = 10, nx = 10,  seed = 11)
  expect_equal(dim(f2), c(10,10))
})
