\name{glm1path}
\alias{glm1path}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fits a path of Generalised Linear Models with LASSO (or L1) penalties, and finds the model that minimises BIC.
}
\description{
Fits a sequence (path) of generalised linear models with LASSO penalties, using an iteratively reweighted local linearisation approach. The whole path of models is returned, as well as the one that minimises BIC. Can handle negative binomial family, even with overdispersion parameter unknown, as well as other GLM families. 
}
 %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
\usage{
glm1path(y, X, family = "negative.binomial", lambdas = NULL,
  penalty = c(0, rep(1, dim(X)[2]-1)), df.max = sum(y > 0), n.lambda = 25, lam.max = NULL,
  lam.min = NULL, k = log(length(y)), b.init = NA, phi.init = NA, phi.iter = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{lambdas}{
An optional vector of LASSO penalty parameters, specifying the path along which models will be fitted. This penalty is applied to parameters as specified in \code{penalty}. By default, a geometric sequence of values will be constructed with \code{n.lambda} values, starting from the intercept model and reducing lambda to 1.e-6 of its original value. Any vector that is provided will be sorted in decreasing order, so that the smallest model (biggest penalty) is fitted first.
}
  \item{penalty}{
A vector to be multiplied by each lambda to make the penalty for each fitted model. The main purpose here is to allow penalties to be applied to some parameters but not others, but it could also be used to change the size of the penalty for some terms as compared to others (e.g. to fit an adaptive LASSO). Must have the same length as the dimension of the model, \code{dim(X)[2]}.
}
  \item{df.max}{
The maximum number of terms that is permitted in the fitted model. Once this threshhold is reached no further fits are attempted. The default break-point is the number of non-zero values in the response vector.
}
  \item{n.lambda}{
The number of models to fit along the path (if not previously specified via \code{lambdas}).
}
  \item{lam.max}{
The maximum value of the LASSO penalty to use along the path of fitted values (if not previously specified via \code{lambdas}).
}
  \item{lam.min}{
The minimum value of the LASSO penalty to use along the path of fitted values (if not previously specified via \code{lambdas}).
}
  \item{k}{
In BIC calculation, this is the value of the penalty per parameter in the fitted model. The default value, \code{log(length(y))}, gives BIC (known to be consistent, for adaptive LASSO), changing it to \code{2} would give AIC (which is not so great in terms of properties).
}
  \item{b.init}{
An initial value for beta for the first model along the fitted path. Default is to fit an intercept model.
}
  \item{phi.init}{
For negative binomial models: An initial value for the overdispersion parameter for the first model along the fitted path. Default is zero (Poisson fit).
}
  \item{phi.iter}{
Number of iterations estimating the negative binomial overdispersion parameter (if applicable) before returning to slope estimation. Default is one step, i.e. iterating between one-step estimates of beta and phi.
}
  \item{\dots}{
Arguments passed to \code{\link{glm1}}.}
}
\details{
This function fits a series of LASSO-penalised generalised linear models, with different values for the LASSO penalty. Largely inspired by the glmnet package. This results in  a path of fitted models, from small ones (with big LASSO penalties) to larger ones (with smaller penalties). Each individual model is fitted using the  \code{\link{glm1}} function, which uses a local linearisation approach as in Osborne et al (2000), nested inside iteratively reweighted (penalised) least squares, and using results from the previous fit as initial estimates. Look it's not the fastest thing going around, try glmnet if you want something faster (and possibly rougher as an approximation). The main advantage of the \code{glm1path} function is that it has been written to accept any glm family argument (although not yet tested beyond discrete data!), and also the negative binomial distribution, which is especially useful for modelling overdispersed counts.

For negative binomial with unknown overdispersion use \code{"negative.binomial"}, or if overdispersion is to be specified, use \code{negative.binomial(theta)} as in the \code{MASS} package. Note that the output refers to phi=1/theta, i.e. the overdispersion is parameterised in output such that the variance is mu+phi*mu^2. Hence values of phi close to zero suggest little overdispersion, values over one suggest a lot.

You can use the \code{residuals} and \code{plot} functions on \code{glm1path} objects in order to compute Dunn-Smyth residuals and a plots of these residuals against linear predictors, as for \code{\link{manyglm}}.
}

\value{An object of class \code{glm1path} with the following components:
\describe{
\item{coefficients}{Vector of model coefficients for the best-fitting model (as judged by BIC)}
\item{lambda}{The value of the LASOS penalty parameter, lambda, for the best-fitting model (as judged by BIC)}
\item{glm1.best}{The glm1 fit for the best-fitting model (as judged by BIC). For what this contains see \code{\link{glm1}}.}
\item{all.coefficients}{A matrix where each column represents the model coefficients for a fit along the path specified by \code{lambdas}.}
\item{lambdas}{A vector specifying the path of values for the LASSO penalty, arranged from largest (strongest penalty, smallest fitted model) to smallest (giving the largest fitted model).}
\item{logL}{A vector of log-likelihood values for each model along the path.}
\item{df}{A vector giving the number of non-zero parameter estimates (a crude measure of degrees of freedom) for each model along the path.}
\item{bics}{A vector of BIC values for each model along the path. Calculated using a penalty on model complexity as specified by input argument \code{k}.}
\item{counter}{A vector counting how many iterations until convergence, for each model along the path.}
\item{check}{A vector of logical values specifying whether or not Karush-Kuhn-Tucker conditions are satisfied at the solution.}
\item{phis}{For negative binomial regression - a vector of overdispersion parameters, for each model along the path.}
\item{y}{The vector of values for the response variable specified as an input argument.}
\item{X}{The design matrix of p explanatory variables specified as an input argument.}
\item{penalty}{The vector to be multiplied by each lambda to make the penalty for each fitted model.}
\item{family}{The family argument specified as input.}
}
}

\references{
Osborne, M.R., Presnell, B. and Turlach, B.A. (2000) On the LASSO and its dual. Journal of Computational and Graphical Statistics, 9, 319-337.
}
\author{
David I. Warton <David.Warton@unsw.edu.au>
}
\seealso{
\code{\link{glm1}}, \code{\link{glm}}, \code{\link{family}}, \code{\link[mvabund:residuals.manyglm]{residuals.manyglm}}, \code{\link[mvabund:plot.manyany]{plot.manyany}}
}
 %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
\examples{
data(spider)
Alopacce <- spider$abund[,1]
X <- cbind(1,spider$x)

# fit a LASSO-penalised negative binomial regression:
ft = glm1path(Alopacce,X)
# have a look at the BICS for all models:
plot(ft$bics~ft$lambdas, log="x")

#the action seems to be at lambda above 0.1, re-do with a minimum lambda at 0.1 and more lambdas:
ft2 = glm1path(Alopacce,X,lam.min=0.1,n.lambda=100)
plot(ft2$bics~ft2$lambdas, log="x")

# return the slope estimates for the best-fitting model:
coef(ft2)

# look at a residual plot:
plot(ft2)

}
\keyword{models}
\keyword{regression}
