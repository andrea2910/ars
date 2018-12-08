library(testthat)
library(ars)

test_that("ARS breaks with incorrect input",{
  expect_error(ARS$new(1)) # first input isn't a function
  expect_error(ARS$new(dnorm, D=c(1,2,3))) # D isn't the correct length
  expect_error(ARS$new(dexp, rate=1, D=c(-5,5))$sample())
  # undefined bound so it should break later on
})

test_that("namedVector breaks with incorrect input",{
  expect_error(namedVector$new(c("a", "b", "c", c(1,2,3))))
  expect_error(namedVector$new(c(12,3), c(2,5,2,4)))
})

test_that("output is correct",{
  expect_length(ARS$new(dnorm, mean=5, sd=2)$sample(n=100), 100)
  expect_length(ARS$new(dexp, D=c(0, 10), rate=2)$sample(n=100), 100)
  expect_length(ARS$new(dunif, D=c(0,1))$sample(n=100), 100)
})

test_that("helper functions are correct",{
  expect_error(calc_deriv(5, 2, 5, 2))
  expect_equal(check_interpolconcave(seq(-3, 3, by=1),
                                     log(dnorm(seq(-3, 3, by=1)))), TRUE)
  expect_type(create_upphull(seq(-3, 3, by=1),
                             log(dnorm(seq(-3, 3, by=1)))), "closure")
  expect_type(create_lowhull(seq(-3, 3, by=1),
                             log(dnorm(seq(-3, 3, by=1)))), "closure")

})

