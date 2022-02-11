context("test-manyglm.R")
source('utils.R')

test_that("bad link function errors", {
  set.seed(100)
  data(spider)
  spiddat <- mvabund(spider$abund)
  #To fit a log-linear model assuming counts are poisson:
  bad_fam <- "avdsfjhsd"
  expect_error(manyglm(spiddat~., data=spider$x, family=bad_fam),"'family'.*")
  expect_error(manyglm(spiddat~., data=spider$x, family=data.frame()),"'family'.*")
})

test_that("spider coefs", {
  data(spider)
  spiddat <- mvabund(spider$abund)
  #To fit a log-linear model assuming counts are poisson:
  glm.spid <- manyglm(spiddat~., data=spider$x, family="poisson")
  expect_is(glm.spid, "manyglm")

  spider_coefs <- coef(glm.spid[1])
  spider_coefsREF <- readRDS("spider_coefs.rds")
  delta <- spider_coefs - spider_coefsREF
  expect_lte(sum(delta^2), 1.e-7)
  spider_summ <- summary(glm.spid,nBoot=1000)
  spider_summREF <- readRDS("poisson_spider_summary.rds")
  deltaWald = spider_summ$coefficients[,1] - spider_summREF$coefficients[,1]
  expect_lte(sum(deltaWald^2), 1.e-7)
  deltaP = spider_summ$coefficients[,2] - spider_summREF$coefficients[,2]
  expect_lte(sum(deltaP^2), 0.003) # note Monte Carlo error possible here
})

test_that('theta estimation method warnings', {
  Y <- sapply(1:4, function(x)  rnbinom(30, 5, 0.4))
  Y <- mvabund(Y)
  # method of moments is not allowed for negative binomial family
  expect_error(manyglm(Y ~ 1, theta.method = 'MM'), 'theta.method.*')
})


test_that("gamma family shape parameter", {
  set.seed(100)
  n <- 1000; shapes <- 1:4; rate <- 1
  Y <- mvabund(sapply(shapes, function(shape) rgamma(n, shape, rate)))
  gamma_glm <- manyglm(Y ~ 1, family="gamma", show.warning = T)
  thetas <- signif(gamma_glm$theta, 5)
  thetasREF <- c(1.0034,2.0542,2.9994,4.0487)
  deltaThetas=thetas-thetasREF
  expect_lte(sum(deltaThetas^2), 1.e-7)
})


