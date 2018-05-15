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

test_that('theta estimation method warnings', {
  Y <- sapply(1:4, function(x)  rnbinom(30, 5, 0.4))
  # method of moments is not allowed for negative binomial family
  expect_error(mvabund(Y ~ 1, theta.method = 'MOMENTS'), 'theta.method.*')
})

gen_gamma <- function(n, shape, rates) {
  sapply(rates, function(rate) rgamma(n, shape, rate))
}

test_that('poisson family', {
  skip('causes seg fault')
  Y <- sapply(1:4, function(x)  rpois(30, x))
  mvabund(Y) -> Y
  manyglm(Y ~ 1 , family = 'poisson') -> m
  print('poisson manyglm summary')
  print(summary(m, show.warning = T))
  junk <- rnorm(1000)
  manyglm(Y ~ junk , family = 'poisson') -> m2
  exp(coef(m2)[1,])
  exp(coef(m2))[2,]
})

test_that("gamma family", {
  skip('gamma family')
  n <- 1000; shape <- 1; rates <- 1:4
  Y <- mvabund(gen_gamma(n, shape, 1:4))
  junk <- rnorm(n)
  # using the log link
  gamma_glm <- manyglm(Y ~ junk, family="gamma", show.warning = T)
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

test_that("gamma family shape parameter", {
  skip('write actual tests')
  n <- 1000; shapes <- 1:4; rate <- 1
  Y <- mvabund(sapply(shapes, function(shape) rgamma(n, shape, rate)))
  print(shapes)
  gamma_glm <- manyglm(Y ~ 1, family="gamma", show.warning = T)
  expect_equal(2*2, 4)
})
test_that("gamma family summary", {
  set.seed(200)
  shape = 1
  n <- 100
  Y <- mvabund(gen_gamma(n, shape, c(1,1,1)))
  # using the log link
  junk <- rnorm(n)
  gamma_glm <- manyglm(Y ~ junk, family="gamma", show.warning = T, k = shape)
  tests <- c('wald','score','LR')
  for (t in tests)
    print(summary(gamma_glm, show.warning = T, test=t))
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

