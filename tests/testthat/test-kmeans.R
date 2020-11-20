




test_that("kmeans works", {

  set.seed(1)
  Y <- matrix(rpois(n = 100 * 30, lambda = 5), nrow = 30, ncol = 100)
  min_mu <- 1e-4
  k <- 3

  size_factors <- matrixStats::colMeans2(Y)
  centers <- Y[, sample(seq_len(ncol(Y)), size = k)] + min_mu

  # Call clustering
  res1 <- run_kmeans(Y, size_factors, centers,
                     min_mu = min_mu, max_iter = 100,
                     tolerance = 1e-8, verbose = FALSE)

  res2 <- run_kmeans_r(Y, size_factors, centers,
                       min_mu = min_mu, max_iter = 100,
                       tolerance = 1e-8, verbose = FALSE)
  expect_equal(res1, res2)
})



test_that("kmeans++ initialization works", {

  set.seed(1)
  Y <- matrix(rpois(n = 100 * 30, lambda = 5), nrow = 30, ncol = 100)
  min_mu <- 1e-4
  k <- 3

  size_factors <- matrixStats::colMeans2(Y)
  centers <- Y[, sample(seq_len(ncol(Y)), size = k)] + min_mu

  set.seed(1)
  res1 <- kmeans_pp_initialization(Y, size_factors, k, min_mu)
  set.seed(1)
  res2 <- kmeans_pp_initialization_r(Y, size_factors, k, min_mu)
  expect_equal(res1, res2)

  head(res1)
  head(res2)
})



