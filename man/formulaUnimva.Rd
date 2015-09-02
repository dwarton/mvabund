\name{formulaUnimva}
\alias{formulaUnimva}
\title{Create a List of Univariate Formulas}
\description{
Create a list of m univariate formulas given a formula with multivariate
response of dimension m.
}
\usage{
formulaUnimva(formula, var.subset, split.x=FALSE, intercept=0,
  allow.noresp=FALSE)
}
\arguments{
\item{formula}{a formula or a mvformula, the elements are not allowed
to be data.frames.}
\item{var.subset}{optional vector of the variable numbers to use.}
\item{split.x}{
logical, whether explanatory terms that are matrices should be split and
each added as a single term. this is useful for plotting formulas.
}
\item{intercept}{
numeric, either \code{1} if an Intercept should be included in the formula 
or \code{0} if there shouldn't be an Intercept in the formula.
}
\item{allow.noresp}{
logical, whether an empty response is allowed (a list with one element
would be returned) or not (would result in an error.)
}
}
% \details{}
\value{
A list containing m formulas with the univariate responses chosen by var.subset.
}
% \references{...}
% \note{...}
\author{
Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{mvformula}}, \code{\link{mvabund}}.
}
\examples{
data(spider)
spiddat <- mvabund(spider$abund)
X <- spider$x

formulaUnimva(spiddat~X)
}
\keyword{models}
\keyword{multivariate}
\keyword{manip}

