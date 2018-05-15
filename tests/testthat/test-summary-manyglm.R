context("test-manyglm.R")
set.seed(1017)
gen_pois <- function(n, m, how_many = 5) {
  replicate(how_many, rpois(n, m))
}

test_that("gamma summary", {
  # generate a response matrix Xij ~ exp(1) distribution
  # ${nVar} cols
  # with one junk predictor from the normal distribution
  shape <- 1; rate <- 1; nVar <- 5
  SEED <- 1001
   <- 'gamma'
  set.seed(SEED)
  Y2 <- mvabund(sapply(rep(1,nVar), function(rate) rgamma(n, shape, rate)))
  junk <- rnorm(n)
  mglmg <- manyglm(Y2 ~ junk, family=FAMILY, show.warning = TRUE, maxiter = 500)

  tests <- c('wald', 'score', 'LR')
  resamp = c('case', 'perm.resid', 'montecarlo', 'pit.trap')
  summaries <- list()
  for (t in tests)
    for (s in resamp) {
      summaries[[paste(FAMILY, t, s)]] <- summary(mglmg,
        test=t,
        resamp = s,
        rep.seed = SEED)
    }
  expect_equal_to_reference(summaries, 'gamma-summaries.rds')
}


