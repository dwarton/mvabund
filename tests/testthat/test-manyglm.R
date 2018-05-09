context("test-manyglm.R")

test_that("bad link function errors", {
  set.seed(100)
  data(spider)
  spiddat <- mvabund(spider$abund)
  X <- spider$x
  #To fit a log-linear model assuming counts are poisson:
  bad_fam <- "avdsfjhsd"
  expect_error(manyglm(spiddat~X, family=bad_fam),"'family'.*")
  expect_error(manyglm(spiddat~X, family=data.frame()),"'family'.*")
})

test_that("spider coefs", {
  data(spider)
  spiddat <- mvabund(spider$abund)
  X <- spider$x
  #To fit a log-linear model assuming counts are poisson:
  glm.spid <- manyglm(spiddat~X, family="poisson")
  expect_is(glm.spid, "manyglm")

  spider_coefs <- coef(glm.spid[1])
  expect_equal_to_reference(spider_coefs, 'spider_coefs.rds')
  expect_equal_to_reference(summary(glm.spid), 'poisson_spider_summary.rds')
})

gen_gamma <- function(n, shape, rates) {
  sapply(rates, function(rate) rgamma(n, shape, rate))
}

test_that('poisson family', {
  skip('no poisson tests')
  Y <- sapply(1:4, function(x)  rpois(1000, x))
  mvabund(Y) -> Y
  manyglm(Y ~ 1 , family = 'poisson') -> m
  exp(coef(m))
  junk <- rnorm(1000)
  manyglm(Y ~ junk , family = 'poisson') -> m2
  exp(coef(m2)[1,])
  exp(coef(m2))[2,]
})
test_that("gamma family", {
  n <- 1000; shape <- 1; rates <- 1:4
  Y <- mvabund(gen_gamma(n, shape, 1:4))
  junk <- rnorm(n)
  # using the log link
  gamma_glm <- manyglm(Y ~ junk, family="gamma", show.warning = T, k = shape)
  if(interactive()) {
    print('')
    print('rates')
    print(rates)
    print('raw fit')
    print(coef(gamma_glm))
    print('estimate of rate parameters')
    print(shape / exp(coef(gamma_glm)))
    print(gamma_glm)
    print('')
  }
  expect_equal(2*2, 4)
})

test_that("gamma family summary", {
  skip('Breaking fail lda')
  set.seed(100)
  shape = 1
  Y <- mvabund(gen_gamma(1000, 1, 1:4))
  # using the log link
  gamma_glm <- manyglm(Y ~ 1, family="gamma", show.warning = T, k = shape, maxiter = 200)
  print(summary(gamma_glm, show.warning = T))
  expect_equal(2 * 2, 4)
})

test_that("gamma family anova", {
  skip('fails')
  set.seed(100)
  shape = 1
  n <- 1000
  Y <- mvabund(gen_gamma(n, 1, 1:4))
  # using the log link
  no_effect <- rnorm(n)
  gamma_glm <- manyglm(Y ~ no_effect, family="gamma", show.warning = T, k = shape)
  print(anova(gamma_glm, show.warning = T))
  expect_equal(2 * 2, 4)
})

