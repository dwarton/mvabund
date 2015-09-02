\name{glm1}
\alias{glm1}
\title{
Fits a Generalised Linear Models with a LASSO (or L1) penalty, given a value of the penalty parameter.
}
\description{
Fits a generalised linear model with a LASSO penalty, using an iteratively reweighted local linearisation approach, given a value of the penalty parameter (lamb). Can handle negative binomial family, even with overdispersion parameter unknown, as well as other GLM families. 
}
 %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
\usage{
glm1(y, X, lambda, family = "negative.binomial", weights = rep(1, length(y)),
     b.init = NA, phi.init = NA, phi.method = "ML", tol = c(1e-08, .Machine$double.eps),
     n.iter = 100, phi.iter = 1)
}
\arguments{
  \item{y}{
A vector of values for the response variable.
}
  \item{X}{
A design matrix of p explanatory variables.
}
  \item{family}{
The family of the response variable, see \code{\link{family}}. Negative binomial with unknown overdispersion can be specified as "negative.binomial", and is the default.
}
  \item{lambda}{
The penalty parameter applied to slope parameters. Different penalties can be specified for different parameters by specifying lamb as a vector, whose length is the number of columns of X. If scalar, this penalty is applied uniformly across all parameters except for the first (assuming that it is an intercept) 
}
  \item{weights}{
Observation weights. These might be useful if you want to fit a Poisson point process model...
}
  \item{b.init}{
Initial slope estimate. Must be a vector of the same length as the number of columns in X.
}
  \item{phi.init}{
Initial estimate of the negative binomial overdispersion parameter. Must be scalar.
}
  \item{phi.method}{
Method of estimating overdispersion.
}
  \item{tol}{
A vector of two values, specifying convergence tolerance, and the value to truncate fitted values at.
}
  \item{n.iter}{
Number of iterations to attempt before bailing.
}
  \item{phi.iter}{
Number of iterations estimating the negative binomial overdispersion parameter (if applicable) before returning to slope estimation. Default is one step, i.e. iterating between one-step estimates of beta and phi.
}
}
\details{
This function fits a generalised linear model with a LASSO penalty, sometimes referred to as an L1 penalty or L1 norm, hence the name glm1. The model is fit using a local linearisation approach as in Osborne et al (2000), nested inside iteratively reweighted (penalised) least squares. Look it's not the fastest thing going around, try \code{glmnet} if you want something faster (and possibly rougher as an approximation). The main advantage of the \code{glm1} function is that it has been written to accept any glm family argument (although not yet tested beyond discrete data!), and also the negative binomial distribution, which is especially useful for modelling overdispersed counts.

For negative binomial with unknown overdispersion use \code{"negative.binomial"}, or if overdispersion is to be specified, use \code{negative.binomial(theta)} as in the \code{MASS} package. Note that the output refers to phi=1/theta, i.e. the overdispersion is parameterised such that the variance is mu+phi*mu^2. Hence values of phi close to zero suggest little overdispersion, values over one suggest a lot.
}

\value{
\item{coefficients}{Vector of parameter estimates}
\item{fitted.values}{Vector of predicted values (on scale of the original response)}
\item{logLs}{Vector of log-likelihoods at each iteration of the model.  The last entry is the log-likelihood for the final fit.}
\item{phis}{Estimated overdispersion parameter at each iteration, for a negative binomial fit.}
\item{phi}{Final estimate of the overdispersion parameter, for a negative binomial fit.}
\item{score}{Vector of score equation values for each parameter in the model.}
\item{counter}{Number of iterations until convergence. Set to Inf for a model that didn't converge.}
\item{check}{Logical for whether the Kuhn-KArush-Tucker conditions are saitsfied.}
}

\references{
Osborne, M.R., Presnell, B. and Turlach, B.A. (2000) On the LASSO and its dual. Journal of Computational and Graphical Statistics, 9, 319-337.
}
\author{
David I. Warton <David.Warton@unsw.edu.au>, Ian W. Renner and Luke Wilson.
}

\seealso{
\code{\link{glm1path}}, \code{\link{glm1}}, \code{\link{glm}}, \code{\link{family}}
}

\examples{
data(spider)
Alopacce <- spider$abund[,1]
X <- cbind(1,spider$x)
#fit a LASSO-penalised negative binomial GLM, with penalty parameter 10:
ft = glm1(Alopacce,X,lambda=10)

plot(ft$logLs) # a plot of the log-likelihood, each iteration to convergence
coef(ft) # coefficients in the final model

}
\keyword{models}
\keyword{regression}
