\name{mvabund}
\alias{mvabund}
\alias{as.mvabund}
\alias{is.mvabund}

\alias{mvabund-class}

\title{Multivariate Abundance Data Objects}
\description{
\code{mvabund} creates an mvabund object.\cr
\code{as.mvabund} attempts to turn its argument into an mvabund object.\cr
\code{is.mvabund} tests if the argument is an mvabund object.\cr
\code{mvabund} is a class of objects for which special-purpose plotting and regression functions have been written in the \code{\link{mvabund-package}}. The above are useful preliminary functions before analysing data using the special-purpose functions. These new functions were written specially for the analysis of multivariate abundance data in ecology, hence the title 'mvabund'.
}
\usage{
mvabund( \dots , row.names=NULL, check.rows=FALSE, check.names=TRUE,
  var.names=NULL, neg=FALSE, na.rm=FALSE )
  
as.mvabund(x)

is.mvabund(x)
}
\arguments{
\item{ \dots }{these arguments are of either the form value or tag = value.
Component names are created based on 
the tag (if present) or the deparsed argument itself.}
\item{row.names}{\code{NULL} or a single integer or character string specifying
a column to be used as row names, or a character
or integer vector giving the row names for the mvabund object.}
\item{check.rows}{
if \code{TRUE} then the rows are checked for consistency of length
and names.
}
\item{check.names}{logical. If \code{TRUE} then the names of the variables are
checked to ensure that they are syntactically valid variable names
and are not duplicated. If necessary they are adjusted (by make.names) so that
they are.}
\item{var.names}{\code{NULL} or a character vector giving the column names for
the mvabund object.}
\item{neg}{character. If \code{FALSE} negative values will cause an error message.}
\item{na.rm}{logical, whether missing values should be removed.}
\item{x}{an \R object.}
}
\details{
It is desirable to convert abundance data to \code{mvabund} objects, to allow
automatic use of all methods for mvabund objects, for example the provided 
methods for plotting, linear and generalised linear model-fitting and inference.  

Some more technical details:  \cr
\code{mvabund} objects always have two dimensions. \cr
\code{mvabund} converts its arguments into an mvabund object. The supplied argument can be a matrix, data frame, a list of vectors, or several vectors as separate arguments. \cr

If elements of the created \code{mvabund} object are not numeric, a warning will be printed. \cr

If \code{row.names} are not supplied, the row names of the \code{mvabund}
object will be \code{NULL}. If the length of row.names does not match the
number of rows or there are duplicates, an error message will result.
} 

\value{
\code{mvabund} and \code{as.mvabund} returns an mvabund object. \cr
\code{is.mvabund} returns \code{TRUE} if \code{x} is a matrix and \code{FALSE}
otherwise.
}
% \note{ }
\author{
Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{unabund}}, \code{\link{mvformula}}, \code{\link{plot.mvabund}},
Also see the \code{\link{mvabund-package}}.
}
\examples{
data(solberg) 

## Create an mvabund object:
solbergdat <- mvabund(solberg$abund)

## Turn solberg$abund into an mvabund object and store as solbergdat:
solbergdat <- as.mvabund(solberg$abund)

## Check if solbergdat  is an mvabund object:
is.mvabund(solbergdat)
}
\keyword{classes}
\keyword{multivariate}

