\name{logLik.manylm}
\alias{logLik.manylm}
\title{Calculate the Log Likelihood}
\description{Calculate the log likelihood of a multivariate linear model.}
\usage{
% logLik.manylm
\method{logLik}{manylm}(object, REML = FALSE, \dots)
}
\arguments{
 \item{object}{a \code{manylm} object from which a log-likelihood value
    should be extracted.}
 \item{\dots}{some methods for this function require additional arguments.}
  \item{REML}{an optional logical value.  If \code{TRUE} the restricted
    log-likelihood is returned, else, if \code{FALSE}, the
    log-likelihood is returned.  Defaults to \code{FALSE}.}
}
\details{
It is assumed that the scale has been estimated
(by maximum likelihood or REML), and all the constants in the
log-likelihood are included.
}
\value{
  Returns an object, say \code{r}, of class \code{logLik} which is a
  number with attributes, \code{attr(r, "df")} (\bold{d}egrees of
  \bold{f}reedom) giving the number of (estimated) parameters in the model.
}
% \references{...}
% \note{...}
% \author{Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.}
% \seealso{\code{\link{...}}.}
\examples{
data(spider)
spiddat <- mvabund(spider$abund)
lm.spider <- manylm(spiddat~., data=spider$x)
logLik(lm.spider)
}
\keyword{models}
\keyword{multivariate}
