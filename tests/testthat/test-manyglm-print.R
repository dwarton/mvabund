options(echo = T)
library(mvabund)
set.seed(100)
sapply(1:3, function(x) rnbinom(100, 3, 0.5)) -> Y
junk1 <- rnorm(100)
junk2 <- rnorm(100)
Y <- mvabund(Y)

mnb <- manyglm(Y ~ junk1) # negative binomial output
sapply(1:10, function(x) rgamma(100, 3, 0.5)) -> Yg
Yg <- mvabund(Yg)
mg <- manyglm(Yg ~ junk1 + junk2,family = 'gamma')
print(mg)

# residual generic function
scatter.smooth(mg$fitted, residuals(mg),
  lpars = list(col = "red"))
