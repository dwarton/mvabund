\name{ridgeParamEst}
\alias{ridgeParamEst}
\title{
Estimation of the ridge parameter
}
\description{
Maximum likelihood estimation of the ridge parameter by cross-validation
}
\usage{
ridgeParamEst(dat, X, weights = rep(1,times=nRows), refs, 
	tol=1.0e-010, only.ridge=FALSE,  doPlot=FALSE,
	col="blue",type="l", \dots) 
}
\arguments{
\item{dat}{the data matrix.}
\item{X}{the design matrix.}
\item{weights}{weights on the cases of the design matrix.}
\item{refs}{
a vector specifying validation group membership. Default is to
construct \code{refs} using a method that is a function of the sample size \code{N}: 
if \code{N<=20}, leave-one-out is used \code{refs=1:N}, if \code{N<=40}, 
10-fold Cross Validation is used where group membership is chosen randomly 
but with equal size groups, otherwise 5-fold CV with random group memberships.
}
\item{tol}{the sensitivity in calculations near zero.}
\item{only.ridge}{logical, whether only the ridge Parameters should be passed back or 
additionally the Cross Validation penalised likelihood.}
\item{doPlot}{logical, whether a plot of -2logL vs a candidate for the ridge parameter should be drawn.}
\item{col}{color of Plot symbols.}
\item{type}{type of Plot symbols.}
\item{\dots}{further plot arguments.}
}
\details{
This function estimates the ridge parameter when applying ridge regularization to a sample correlation matrix of residuals. The ridge parameter is estimated to maximize the normal likelihood as estimated via cross validation (Warton 2008). 
}
\value{
A list with the following component:
\item{ridgeParameter}{
the estimated ridge parameter
}
If \code{only.ridge=FALSE} the returned list additionally contains the element:
\item{minLL}{
the minimum of the negative log-likelihood}.
}
\references{
Warton D.I. (2008). Penalized normal likelihood and ridge regularization of
correlation and covariance matrices. \emph{Journal of the American
Statistical Association} 103, 340-349.
}
\author{
David Warton <David.Warton@unsw.edu.au> and Ulrike Naumann.
}
\seealso{
\code{\link{manylm}}
}
\examples{
data(spider)
spiddat <- mvabund(spider$abund)
X <- as.matrix(spider$x)

ridgeParamEst(dat = spiddat, X = model.matrix(spiddat~X))
}
\keyword{models}
