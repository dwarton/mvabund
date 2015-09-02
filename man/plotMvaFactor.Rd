\name{plotMvaFactor}
\alias{plotMvaFactor}
\title{Draw a Mvabund Object split into groups. }
\description{Draw the \code{mvabund} object \code{x} but split the data into
groups according to the grouping variable \code{y}.}
\usage{
plotMvaFactor(x, y, type="p", main="Abundance", n.vars= min(12,NCOL(x)), 
  transformation="log", legend=TRUE, ...)
}
\arguments{
\item{x}{a \code{mvabund object}, a matrix with multivariate abundance data.}
\item{y}{a factor or a data.frame with factors, non-factor columns in a
data.frame are ignored.}
\item{type}{what type of plot should be drawn, allowed types are "p" for
scatterplot, "bx" for boxplot and "n" for no plot. Other types, as used
in \code{par} are NOT allowed.}
\item{main}{the title of the plot, see \code{plot}.}
\item{n.vars}{the number of variables to include in the plot.}
\item{transformation}{an optional transformation, "no" = untransformed,
"sqrt"=square root transformed, "log" (default)=log(Y/min+1) transformed,
"sqrt4" =4th root transformed.}
\item{legend}{logical, whether a legend should be added to the plot.}
\item{\dots}{arguments to be passed to or from other methods.}
}
\details{
For each variable in y that is a factor, a plot is drawn. When boxplots are drawn
the colors, that can be supplied by \code{col} are used to display different
factor levels.
For scatterplots it is also possible to use the plotting symbols, specified by
\code{pch} for that. \cr

If the colors and for scatterplots the plotting symbols are not supplied,
they will be automatically generated. However, the plotting symbols will only
be automatically used in this way if there are up to seven different levels.

If colors or the plotting symbols are supplied, but the number of factor levels
is bigger than the the number of different values, they will be replicated.

Sometimes the legends might be only partially visible, especially when the width
of the graphics device is too small. To fix this, create a graphics device with
a larger width (see help("device") for on available devices and their details)
and then repeat the
\code{plotMvaFactor} command.
}
% \value{...}
\references{
Warton, D. I. ( )
\emph{Raw data graphing: an informative but under-utilised tool 
for the analysis of multivariate abundances}, , .
}
% \note{...}
\author{
Ulrike Naumann, Yi Wang, Stephen Wright and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{plot.mvabund}}.
}
\examples{
require(graphics)

## Plot an Environment Factor vs Abundance plot
data(spider)
spiddat <- mvabund(spider$abund)
X <- spider$x

## Create a Environmental factor where TRUE=Sand, FALSE=No Sand)
X <- as.factor(X[,2]>0)
plotMvaFactor(x=spiddat, y=X)
}
\keyword{hplot}

