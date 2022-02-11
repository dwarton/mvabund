\name{best.r.sq}
\alias{best.r.sq}
\title{
Use R^2 to find the variables that best explain a multivariate response.
}
\description{
Finds the subset of explanatory variables in a formula that best explain the variation in a multivariate response, as measured by a chosen definition of R^2. Modifications are included for high dimensional data, such as multivariate abundance data in ecology. 
}
\usage{
best.r.sq(formula, data = parent.frame(), subset, var.subset,
  n.xvars= min(3, length(xn)), R2="h", \dots)
}
\arguments{
\item{formula}{a mvformula, a multivariate formula.}
\item{data}{optional, the data.frame (or list) from which the variables in
formula should be taken.}
\item{subset}{an optional vector specifying a subset of observations to be
used in the fitting process.}
\item{var.subset}{an optional vector specifying the subset of the responses
to be used.}
\item{n.xvars}{the number of independent variables with the highest average
R^2 that should be found.
}
\item{R2}{
 the type of R^2 (correlation coefficient) that should be shown, possible values are: \cr
 "h" = Hooper's R^2 = tr(SST^(-1)SSR))/p   \cr
% "t" = using the trace R^2 =tr(SSR)/tr(SST) \cr
 "v" = vector R^2 = det(SSR)/det(SST)\cr
 "n" = none
Note that for a univariate response, all of these are equivalent to the ordinary 
product-moment correlation coefficient. }
\item{\dots}{further arguments that are passed on to lm.} }

\details{ 
\code{best.r.sq} finds the \code{n.xvars} influence variables obtained by a 
forward selection in a multivariate linear model given by \code{formula}. \cr 
Only the response variables given by \code{var.subset} are considered. However, 
if \code{var.subset} is \code{NULL} all response variables are considered.\cr 
Interactions are excluded from the search mechanism, however the indices that are 
returned correspond to the indices in the model. This function is intended as an 
exploratory tool which can be used for example in plotting, and is not intended 
as a tool for formal model selection. % add an option 'method' that allows to 
choose 'all possible subsets' % instead of 'forward selection' (which is used at 
the moment) } 

\value{This function returns a list consisting of:
\describe{
\item{xs}{a vector of indices of independent variables with the greatest
explanatory power, as previously.}
\item{r2Step}{a vector of total R^2 from sequential model fits including each of
the model terms identified in \code{xs}.}
\item{r2Matrix}{a matrix containing the total R^2 for each term in the model at
each addition step (steps in columns and model terms in rows).}
}
} 

\author{
Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
% \seealso{}
\examples{
data(spider)
spiddat <- mvabund(spider$abund)
X <- as.matrix(spider$x)

best.r.sq( spiddat~X )
}
\keyword{htest}
\keyword{models}
\keyword{regression}
\keyword{multivariate}
