\name{extend.x.formula}
\alias{extend.x.formula}
\title{Extend a Formula to all of it's Terms}
\description{extend a compact formula to all of it's terms as they are interpreted}
\usage{
extend.x.formula(formula, extend.term=TRUE, return.interaction=TRUE)
}
\arguments{
\item{formula}{a model formula.}
\item{extend.term}{logical. If \code{TRUE} terms that refer to multiple variables
 are split into it's multiple terms.}
\item{return.interaction}{logical. Whether a list containing the new formula and 
a vector containing logical values with information about interactions should be returned
or only the new formula.}
}
% \details{...}
\value{
If \code{return.interaction} is \code{TRUE} a list containing the components:
\item{formula}{the new formula.}
\item{is.interaction}{logical, vector giving information 
whether the corresponding formula term is an interaction or not.}
}
% \references{...}
% \note{...}
\author{Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{mvformula}}, \code{\link{formulaUnimva}}, \code{\link{plot.mvformula}},
\code{\link{best.r.sq}}.
}
\examples{
data(spider)
spiddat <- mvabund(spider$abund)
X <- spider$x

foo <- mvformula(spiddat~ X[,1]*X[,2]+log(X[,3]))
extend.x.formula(foo)
}
\keyword{models}
