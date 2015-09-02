\name{mvformula}
\alias{mvformula}
\alias{as.mvformula}
\alias{is.mvformula}

\alias{mvformula-class}

\title{ Model Formulae for Multivariate Abundance Data }
\description{
\code{mvformula} is a method to create an object of class \code{mvformula} \cr
\code{as.mvformula} is a function to turn a formula into a \code{mvformula} \cr
\code{is.mvformula} tests if its argument is a \code{mvformula} object.
\code{mvformula} is a class of objects for which special-purpose plotting and regression functions have been written in the \code{\link{mvabund-package}}. The above are useful preliminary functions before analysing data using the special-purpose functions. These new functions were written specially for the analysis of multivariate abundance data in ecology, hence the title 'mvabund'.
}
\usage{
mvformula(x)
% \usage{formulaMva(x)}
as.mvformula(x)
is.mvformula(x)
}
\arguments{
\item{x}{an \R object.}
}
\details{
\code{mvformula} is a method to create an object of class \code{mvformula}
If the response of the resulting formula is not a \code{mvabund} object, a warning
will be printed. 
\code{as.mvformula} is a function to turn a formula into a \code{mvformula}
if the response in x is a data.frame or an unsuitable object the conversion 
will fail.
}
\value{a formula mvabund for \code{mvformula} and \code{as.mvformula} 
a logical value indicating whether \code{x} is a \code{mvformula} object }
\author{
Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{mvabund}}.
\code{\link{manylm}}.
}
\examples{
require(graphics)

data(spider)
spiddat <- mvabund(spider$abund)
X <- spider$x
 
## Create a formula for multivariate abundance data:
foo <- mvformula( spiddat~X )

## Check whether foo is a mvformula:
is.mvformula(foo)

## Plot a mvformula:
plot(foo)

}
\keyword{models}
\keyword{multivariate}
