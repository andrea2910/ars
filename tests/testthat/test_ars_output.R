library(testthat)
library(ars)

context("ARS behaves correctly")

test_that("ARS breaks with incorrect input",{
  expect_error(ars$new(1)) # first input isn't a function
  expect_error(ars$new(dnorm, D=c(TRUE,FALSE)))
})

test_that("output is correct",{
  expect_length(ars$new(dnorm, mean=5, sd=2)$sample(n=100), 100)
  expect_length(ars$new(dexp, D=c(0, 10), rate=2)$sample(n=100), 100)
  expect_length(ars$new(dunif, D=c(0,1))$sample(n=100), 100)
})

test_that("ARS breaks with not log-concave functions or incorrect bounds",{
  expect_error(ars$new(dexp, rate=5, D=c(-5,5))$sample()) #also has a warning
  expect_error(ars$new(dunif, D=c(-1,-0.5))$sample())
  expect_error(ars$new(dbeta, D=c(0,1), shape1=0.5, shape2=0.5)$sample())

  f <- function(x) dnorm(x, mean=5, sd=2) - exp(-5)
  expect_error(ars$new(f)$sample(n=100))

  f2 <- function(x) dbeta(x, shape1=0.25, shape2=10)
  expect_error(ars$new(f2, D=c(0,1))$sample(n=100))
})
