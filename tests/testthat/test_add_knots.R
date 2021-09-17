test_that("testing the add_knots", {
  # test vectore data class
  f_vec <- rbetafda(n = 1, nx = 30, seed = 3)
  expect_equal(class(add_knots(f = f_vec, knots = c(0,length(f_vec)), L = 2)), "list")
  KS_1 <- add_knots(f = f_vec, knots = c(0,length(f_vec)), L = 2)
  expect_equal(class(add_knots(f = f_vec, knots = KS_1[[1]], L = 2)), "list")
  # test matrix data class
  f_mat <- rbetafda(n = 3, nx = 10, seed = 21)
  expect_warning(add_knots(f = f_mat, knots = c(0, dim(f_mat)[2]), L = 5), "There are only 2 knots. Reduce L or M.")
})


# test_that("testing warning",{
#   f_mat <- rbetafda(n = 3, nx = 10, seed = 11)
#   expect_warning(add_knots(f = f_mat, knots = c(0, dim(f_mat)[2]), L = 5), "There are only 2 knots. Reduce L or M.")
# })
