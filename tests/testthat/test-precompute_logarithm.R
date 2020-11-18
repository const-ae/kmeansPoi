
test_that("conversion back and forth works", {

  d <- pi
  e <- get_exponent_as_integer_from_double(d)
  s <- get_significand_as_integer_from_double(d)
  expect_equal(d, fuse_ints_for_double(e, s), tolerance = 0.01)

  e <- 17
  s <- 11
  d <- fuse_ints_for_double(e, s)
  expect_equal(e, get_exponent_as_integer_from_double(d))
  expect_equal(s, get_significand_as_integer_from_double(d))

})


test_that("convertion")




test_that("log_approx works approximately", {
  log(1e300)
  log_approx(1e300)


  x <- 2^runif(n = 1000, -50, 50)
  plot(log(x), log_approx2(x)); abline(0,1)

  # Constant absolute error
  plot(x, log(x) - log_approx2(x), log ="x")
  # Not constant relative error
  plot(x, log(x) / log_approx2(x), log ="xy", ylim = c(0.99, 1.01))

  expect_equal(log(x), log_approx(x), tolerance = 0.01)
})

x <- 2^runif(n = 1e6, -50, 50)
bench::mark(benchmark_log(x),
            benchmark_log_approx(x),
            benchmark_log_approx2(x),
            check = FALSE)




