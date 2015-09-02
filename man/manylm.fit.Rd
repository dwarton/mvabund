\name{manylm.fit}
\alias{manylm.fit}
\alias{manylm.wfit}

\title{workhose functions for fitting multivariate linear models}
\description{
  These are the workhorse functions called by \code{\link{manylm}} used
  to fit multivariate linear models.  These should usually \emph{not} be used
  directly unless by experienced users.
}
\usage{
manylm.fit(x, y, offset = NULL, tol=1.0e-010, singular.ok = TRUE, \dots)
manylm.wfit(x, y, w, offset = NULL, tol=1.0e-010, singular.ok = TRUE, \dots)
}
\arguments{
  \item{x}{design matrix of dimension \code{n * p}.}
  \item{y}{matrix or an \code{mvabund} object of observations of dimension \code{n*q}.}
  \item{w}{vector of weights (length \code{n}) to be used in the fitting
    process for the \code{manylm.wfit} functions.  Weighted least squares is
    used with weights \code{w}, i.e., \code{sum(w * e^2)} is minimized.}
  \item{offset}{numeric of length \code{n}).  This can be used to
    specify an \emph{a priori} known component to be included in the
    linear predictor during fitting.}
  \item{tol}{tolerance for the \code{qr} decomposition.  Default
    is 1.0e-050.}
  \item{singular.ok}{logical. If \code{FALSE}, a singular model is an
    error.}
  \item{\dots}{currently disregarded.}
}
% \details{}
\value{
  a list with components
  \item{coefficients}{\code{p} vector}
  \item{residuals}{\code{n} vector or matrix}
  \item{fitted.values}{\code{n} vector or matrix}
  % \item{effects}{(not null fits)\code{n} vector of orthogonal single-df
  %   effects.  The first \code{rank} of them correspond to non-aliased
  %   coeffcients, and are named accordingly.}
  \item{weights}{\code{n} vector --- \emph{only} for the \code{*wfit*}
    functions.}
  \item{rank}{integer, giving the rank}
  \item{qr}{(not null fits) the QR decomposition.}
  \item{df.residual}{degrees of freedom of residuals}
  \item{hat.X}{the hat matrix.}
  \item{txX}{the matrix \code{(t(x)\%*\%x)}.}
}
% \references{...}
% \note{...}
\author{
Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
 \code{\link{manylm}} 
}
% \examples{}
\keyword{regression}
\keyword{multivariate}
\keyword{array}

