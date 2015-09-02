\name{shiftpoints}
\alias{shiftpoints}
\title{Calculate a shift for plotting overlapping points}
\description{Calculate a shift to add to overlapping points in plots for better visibility}
\usage{
shiftpoints(x, y, sh=( max(x)-min(x))/100, centered=TRUE, method=1, reg=6,
na.rm=TRUE)
}
\arguments{
\item{x}{a data matrix or numeric vector for use in a plot.}
\item{y}{a data matrix or numeric vector for use in a plot.}
\item{sh}{the maximum total shift.}
\item{centered}{logical, whether the shift is centered at 0, if \code{FALSE}
the shift will be positive only.}
\item{method}{numeric, can have the value 1 or 2, see Details.}
\item{reg}{numeric, see Details.}
\item{na.rm}{logical, indicating whether missing values should be removed.}
}
\details{
This function is similar to \code{\link{jitter}} but is defines for points in
a two-dimensional plot. In contrast to \code{jitter} only the points with ties
have a shift different from 0. The method to calculate the shift is therefore 
not based on random numbers. \cr
If \code{method=1} (default) the individual shift will be selected so that the
shift range is \code{sh}, without regard of the number of overlapping points \cr
\code{method=2} means that for up to \code{reg} overlapping values a fixed
shift of sh/reg is used
}
\value{
Returns an array of shift values with the same dimension
as \code{x}. 
}
% \references{...}
% \note{...}
\author{
Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{ \code{\link{plot.mvabund} }, \code{\link{plot.mvformula}},
\code{\link{jitter}}. }
\examples{
shiftpoints( x=c(1:5, 1:10), y=c(2:6, 1:10) )
}
\keyword{dplot}

