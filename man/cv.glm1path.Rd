\name{cv.glm1path}
\alias{cv.glm1path}
\title{
Fits a path of Generalised Linear Models with LASSO (or L1) penalties, and finds the best model by corss-validation.
}
\description{
Fits a sequence (path) of generalised linear models with LASSO penalties, using an iteratively reweighted local linearisation approach. The whole path of models is returned, as well as the one that minimises predictive log-likelihood on random test observations. Can handle negative binomial family, even with overdispersion parameter unknown, as well as other GLM families. 
}
 %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
\usage{
cv.glm1path(object, block = NULL, best="min", plot=TRUE, prop.test=0.2, n.split = 10,
    seed=NULL, show.progress=FALSE, \dots)
}
\arguments{
  \item{object}{
Output from a \code{\link{glm1path}} fit.
}
  \item{block}{
A factor specifying a blocking variable, where training/test splits randomly assign blocks of observations to different groups rather than breaking up observations within blocks. Default (\code{NULL}) will randomly split rows into test and training groups.
}
  \item{best}{
How should the best-fitting model be determined? \code{"1se"} uses the one standard error rule, \code{"min"} (or any other value) will return the model with best predictive performance. WARNING: David needs to check se calculatios...
}
  \item{plot}{
Logical value indicating whether to plot the predictive log-likelihood as a function of model complexity.
}
  \item{prop.test}{
The proportion of observations (or blocks) to assign as test observations. Default value of 0.2 gives a 80:20 training:test split.
}
  \item{n.split}{
The number of random training/test splits to use. Default is 10 but the more the merrier (and the slower).
}
  \item{seed}{
A vector of seeds to use for the random test/training splits. This is useful if you want to be able to exactly replicate analyses, without Monte Carlo variation in the splits. Default will not used fixed seeds.
}
  \item{show.progress}{Logical argument, if TRUE, console will report when a run for a seed has been completed. This option has been included because this function can take yonks to run on large datasets.
  }
  \item{\dots}{Further arguments passed through to \code{\link{glm1path}}.
  }
}

\details{
This function fits a series of LASSO-penalised generalised linear models, with different values for the LASSO penalty, as for \code{\link{glm1path}}. The main difference is that the best fitting model is selected by cross-validation, using \code{n.test} different random training/test splits to estimate predictive performance on new (test) data. Mean predictive log-likelihood (per test observation) is used as the criterion for choosing the best model, which has connections with the Kullback-Leibler distance. The \code{best} argument controls whether to select the model that maximises predictive log-likelihood, or the smallest model within 1se of the maximum (the '1 standard error rule').

All other details of this function are as for \code{\link{glm1path}}.
}
\value{
\item{coefficients}{Vector of model coefficients for the best-fitting model (as judged by predictive log-likelihood)}
\item{lambda}{The value of the LASOS penalty parameter, lambda, for the best-fitting model (as judged by predictive log-likelihood)}
\item{glm1.best}{The glm1 fit for the best-fitting model (as judged by predictive log-likelihood). For what this contains see \code{\link{glm1}}.}
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
\item{ll.cv}{The mean predictive log-likelihood, averaged over all observations and then over all training/test splits.}
\item{se}{Estimated standard error of the mean predictive log-likelihood.}
}

\references{
Osborne, M.R., Presnell, B. and Turlach, B.A. (2000) On the LASSO and its dual. Journal of Computational and Graphical Statistics, 9, 319-337.
}
\author{
David I. Warton <David.Warton@unsw.edu.au>
}
\seealso{
\code{\link{glm1path}}, \\code{\link{glm1}}, \code{\link{glm}}, \code{\link{family}}
}
 %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
\examples{
data(spider)
Alopacce <- spider$abund[,1]
X <- cbind(1,spider$x)

# fit a LASSO-penalised negative binomial regression:
ft = glm1path(Alopacce,X,lam.min=0.1)
coef(ft)

# now estimate the best-fitting model by cross-validation:
cvft = cv.glm1path(ft)
coef(cvft)

}
\keyword{models}
\keyword{regression}
