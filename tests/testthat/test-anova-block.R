context("test-anova-block.R")
test_that("manyglm block anova ", {
  # from the example(anova.manyglm) 
  ### Load the Tasmania data set
  rep.seed <- TRUE
  data(Tasmania)

  ## Visualise the effect of treatment on copepod abundance
  tasm.cop <- mvabund(Tasmania$copepods)
  treatment <- Tasmania$treatment
  block <- Tasmania$block

  ## Fitting predictive models using a negative binomial model for counts:
  tasm.cop.nb <- manyglm(tasm.cop ~ block*treatment, family="negative.binomial")
  # test just using the normal distribution
  # note this would be inappropriate to do normally but we are trying to test block resampling
  tasm.cop <- manylm(tasm.cop ~ block*treatment)

  ## Testing hypotheses about the treatment effect and treatment-by-block interactions,
  ## using a Wald statistic and 199 resamples (better to ramp up to 999 for a paper):
  aglm <- anova(tasm.cop.nb, nBoot=199, test="LR", show.time='none', block = block, rep.seed=T)
  alm <- anova(tasm.cop, nBoot=199, test="LR", block = block, rep.seed=T)
  expect_equal_to_reference(cbind(aglm$table[,3],alm$table[,3]), 'block_resampling.rds')
})
