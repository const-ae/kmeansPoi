
test_that("simple tests", {
  expect_equal(poisson_deviance(3, 2.2),  poisson_deviance_r(3, 2.2))
  expect_equal(poisson_deviance(3, 100),  poisson_deviance_r(3, 100))
  expect_equal(poisson_deviance(15L, 0.1),  poisson_deviance_r(15L, 0.1))
  y <- rpois(n = 100, lambda = 14.1)
  mu <- exp(rnorm(n = 100, mean = 0.2))
  expect_equal(poisson_deviance(y, mu),  poisson_deviance_r(y, mu))
  expect_equal(poisson_deviance(1, 0.99999999999994),  poisson_deviance_r(1, 0.99999999999994))

})

poisson_deviance(3, 2.2)
poisson_deviance_opt(3, 2.2)
poisson_deviance_r(3, 2.2)
