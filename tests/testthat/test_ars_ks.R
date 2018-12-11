library(testthat)
library(ars)

set.seed(5)
context("ARS passes Kolmogorov-Smirnov Tests for various dist.")

  custom_check <- function(sample1, sample2, alpha=0.05){
    pval <- ks.test(sample1, sample2, alternative="two.sided")$p.value
    if(pval <= alpha) return(FALSE)
    else return(TRUE)
  }

test_that("ARS passes KS test once", {

  n = 5000

  norm_ars <- ARS$new(dnorm)$sample(n)
  norm_r <- rnorm(n)

  unif_ars <- ARS$new(dunif, D=c(0,1))$sample(n)
  unif_r <- runif(n)

  exp_ars <- ARS$new(dexp, D=c(0, 10), rate=2)$sample(n=n)
  exp_r <- rexp(n, rate=2)

  norm_ars_sd50 <- ARS$new(dnorm, mean=1000, sd=50)$sample(n)
  norm_r_sd50 <- rnorm(n, mean=1000, sd=50)

  norm_ars_sd10 <- ARS$new(dnorm, mean=10, sd=10)$sample(n)
  norm_r_sd10 <- rnorm(n, mean=10, sd=10)

  expect_equal(custom_check(norm_ars, norm_r), TRUE)
  expect_equal(custom_check(norm_ars_sd50, norm_r_sd50), TRUE)
  expect_equal(custom_check(norm_ars_sd10, norm_r_sd10), TRUE)
  expect_equal(custom_check(unif_ars, unif_r), TRUE)
  expect_equal(custom_check(exp_ars, exp_r), TRUE)

})

# code below takes too long to run, but idea is the same as above.
# test_that("ARS passes KS test out 90% of 100 times. n=1000", {
#
#   n = 1000
#   pass_norm <- sum(sapply(1:100, function(x) custom_check(ARS$new(dnorm)$sample(n),
#                                                       rnorm(n))))
#
#   pass_unif <- sum(sapply(1:100, function(x) custom_check(ARS$new(dunif, D=c(0,1))$sample(n),
#                                                           runif(n))))
#
#   pass_exp <- sum(sapply(1:100, function(x) custom_check(ARS$new(dexp, D=c(0, 10), rate=2)$sample(n),
#                                                          rexp(n, rate=2))))
#
#   pass_norm_m1000_sd50 <- sum(sapply(1:100, function(x) custom_check(ARS$new(dnorm, mean=1000, sd=50)$sample(n),
#                                                                      rnorm(n, mean=1000, sd=50))))
#
#   pass_norm_m10_sd10 <- sum(sapply(1:100, function(x) custom_check(ARS$new(dnorm, mean=10, sd=10)$sample(n),
#                                                                    rnorm(n, mean=10, sd=10))))
#
#   #pass_norm_m3000_sd1 <- sum(sapply(1:100, function(x) custom_check(ARS$new(dnorm, mean=3000, sd=1)$sample(n),
#    #                                                                rnorm(n, mean=3000, sd=1))))
#
#   expect_gte(pass_norm, 90)
#   expect_gte(pass_unif, 90)
#   expect_gte(pass_exp, 90)
#   expect_gte(pass_norm_m1000_sd50, 90)
#   expect_gte(pass_norm_m10_sd10, 90)
# })
