\name{deviance.manylm}
\alias{deviance.manylm}
% for manyglm the R function deviance.glm works
\title{Model Deviance}
\description{Returns the deviance of a fitted multivariate model object
for abundance data.}
\usage{
% deviance.manylm
\method{deviance}{manylm}(object, na.action="na.omit", \dots)
}
\arguments{
\item{object}{the manylm object.}
\item{na.action}{
how to deal with \code{NA} values. Can be one of 
"na.omit", "na.exclude", "na.fail", \code{NULL}
}
\item{\dots}{additional optional arguments.}
}
% \details{ To extract deviances for fitted models.}
\value{The value of the deviance extracted from the object \code{object}.}
% \references{...}
% \note{...}
% \author{Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.}
\seealso{
\code{\link{manylm}}.
}
\examples{
data(spider)
spiddat <- mvabund(spider$abund)

## Calculate the deviance:
deviance(manylm(spiddat~., data=spider$x))
}
\keyword{models}
\keyword{multivariate}
