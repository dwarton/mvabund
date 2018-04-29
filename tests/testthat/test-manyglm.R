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
test_that("scrap", {
  #expect_equal(2 * 2, 4)
})

