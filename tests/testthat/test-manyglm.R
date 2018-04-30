context("test-manyglm.R")

test_that("link function erros", {
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
})


