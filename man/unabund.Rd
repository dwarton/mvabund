\name{unabund}
\alias{unabund}
\title{ Remove the mvabund Class Attribute }
\description{
Change an mvabund object to a non-mvabund object.
}
\usage{
unabund(x) 
}
\arguments{
\item{x}{an mvabund object that should be transformed into a matrix.}
}
\details{\code{\link{unabund}} doesn't convert \code{x} but only removes the
mvabund class attribute.
}
\value{A matrix if \code{x} is an \code{mvabund} object otherwise \code{x} .}

\author{
Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{mvabund}}.
\code{\link{as.mvabund}}.
\code{\link{is.mvabund}}.
}
\examples{
## Create an mvabund object:
abundances <- as.mvabund(matrix(1:20,5,4))

## Restore the original object:
mat <- unabund(x=abundances)
mat
}
\keyword{sysdata}
\keyword{manip}
\keyword{classes}
