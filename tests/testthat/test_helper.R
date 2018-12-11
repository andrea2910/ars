library(testthat)
library(ars)

context("helper functions")

test_that("helper functions are correct",{
  expect_error(calc_deriv(5, 2, 5, 2))
  expect_equal(check_interpolconcave(seq(-3, 3, by=1),
                                     log(dnorm(seq(-3, 3, by=1)))), TRUE)
  expect_type(create_upphull(seq(-3, 3, by=1),
                             log(dnorm(seq(-3, 3, by=1)))), "closure")
  expect_type(create_lowhull(seq(-3, 3, by=1),
                             log(dnorm(seq(-3, 3, by=1)))), "closure")

})


test_that("namedVector breaks with incorrect input",{
  expect_error(namedVector$new(c("a", "b", "c", c(1,2,3))))
  expect_error(namedVector$new(c(12,3), c(2,5,2,4)))
})
