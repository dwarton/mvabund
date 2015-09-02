\name{boxplot.mvabund}
\alias{boxplot.mvabund}
\alias{boxplot.mvformula}
\title{Boxplots for multivariate abundance Data}
\description{
Draw Boxplots of \code{mvabund} or \code{mvformula} Objects
}
\usage{
\method{boxplot}{mvabund}(x, y, range=1.5, names=NULL, at=NULL, 
                 n.vars=min(12,NCOL(x)), overall.main="Boxplot", 
                 var.subset=NA, transformation="log", ...)

\method{boxplot}{mvformula}(
    x, n.vars=12, overall.main="", var.subset=NA, ...)
}
\arguments{
\item{x}{
for the \code{mvabund} method \code{x} specifies the data from
which the boxplots are to be produced. 
This can be either a numeric vector, or a single list containing such vectors. 
Additional unnamed arguments specify further data as separate vectors 
(each corresponding to a component boxplot). 
NAs are allowed in the data. \cr
For the default method, unnamed arguments are additional data vectors 
(unless \code{x} is a list when they are ignored), and named arguments are
arguments and graphical parameters to be passed to in addition to the ones given
by argument pars (and override those in pars). \cr
For the \code{mvformula} method, a formula, such as \code{y ~ grp}, where
y is a numeric mvabund object of data values to be split into groups according
to the grouping variable grp (a factor). 
}
\item{y}{for the \code{mvabund} method \code{y} can be an additional \code{mvabund} object, if \code{x} isn't a list.}
\item{range}{this determines how far the plot whiskers extend out from the box. If range is positive, the whiskers extend to the most extreme data point which is no more than range times the interquartile range from the box. A value of zero causes the whiskers to extend to the data extremes.}
\item{names}{only available for the \code{mvabund} method: group labels which will be printed under each boxplot.}
\item{at}{only available for the \code{mvabund} method: numeric vector giving the locations where the boxplots should be drawn; defaults to \code{1:n} where \code{n} is the number of boxes.}
\item{n.vars}{the number of variables to include in the plot.}
\item{overall.main}{a character to display as title for every window.}
\item{var.subset}{a numeric vector of indices indicating which variables of the mvabund.object should be included on the plot.}
\item{transformation}{ an optional transformation, (ONLY) for the \code{mvabund} method. Note, that for the \code{mvabund} method \code{transformation} must be used 
instead of \code{log}.\cr
Available values are:\cr
"no" = untransformed, "sqrt"=square root transformed,
"log" (default)=log(Y/min+1) transformed, "sqrt4" =4th root transformed.
}
\item{\dots}{
for the \code{mvformula} method, named arguments to be passed to the \code{plot.mvformula} method. Some arguments that are available for the \code{mvabund} method, are not available in \code{plot.mvformula} and can therefore not available in the \code{mvformula} method.

For the \code{mvabund} method, unamed arguments are additional data of vectors or matrices or \code{mvabund} objects, (unless \code{x} is a list when they are ignored),and named arguments are arguments and graphical parameters to be passed in addition to the ones given by argument \code{pars} (and override those in \code{pars}).
}
}
\details{ 
The function \code{boxplot.mvabund} allows simultaneous construction of many
variables on a single figure. Thus a good comparative overview about the
distribution of abundances for several species can be obtained. \cr
There are several ways in which this function can be used. 
If one \code{mvabund} object, either named \code{x} or \code{y} 
or not names, is passed, it will be drawn on one plot and abundances can be
compared over several variables.\cr
If two \code{mvabund} objects, named \code{x} and \code{y} are 
passed for plotting, they will be shown on
one plot, showing for each species the abundances of both objects directly
one below the other. \cr
If more than two \code{mvabund} objects are passed, each of them will be
plotted separately. \cr
Additionally, it is possible to specify \code{x} as a list of \code{mvabund} objects.
Each of them will be plotted separately and any further \code{mvabund} data will 
be ignored, regardless if it is specified as \code{y} or unnamed.

The function \code{boxplot.mvformula} can be used to draw boxplots of a \code{mvabund}
object in dependence of explanatory variables. The explanatory variables can be both 
numerical values as well as factor variables. If the formula contains both of them,
there will be separate plots for the terms with numerical values and the terms
with factor variables, displayed on separate windows.

The arguments \code{plot}, \code{varwidth} and \code{add}, which are availabe in the default method of \code{boxplot}, are not available for the \code{mvabund} and \code{mvformula} methods. The argument \code{horizontal} is not available for the \code{mvabund} method. \cr
A number of other arguments like \code{at} and \code{names} are only available for the 
\code{mvabund} method.
}
\section{Warning}{
The argument \code{log}, that is available in most plotting functions can not be used 
for plotting \code{mvabund} or \code{mvformula} objects. Instead use \code{transformation} for the \code{mvabund} method and for the \code{mvformula} method
include transformations in the formula.
}

\value{In contrast to the default method (boxplot.default) nothing will be returned.
These functions are only used for drawing the plots.
}
\references{
Warton, D. I. ( )
\emph{Raw data graphing: an informative but under-utilised tool
for the analysis of multivariate abundances}, , .
}
% \note{  }
\author{
Ulrike Naumann, Yi Wang, Stephen Wright and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{plot.mvabund}}.
}
\examples{
require(graphics)

#### Basic Use ####
data(spider)
spiddat <- spider$abund
X <- spider$x

## Create the mvabund object:
spiddat <- mvabund(spiddat)

## Draw a boxplot for a mvabund object:
boxplot(spiddat)

## the same plot could be done by
plot(spiddat,type="bx")


#### Advanced Use ####
data(solberg)
solbdat <- mvabund(solberg$abund)
treatment<- solberg$x

# create pch type and colour vectors
treat.pch <- treat.col <- unclass(treatment)

# Boxplot for data
plot.mvabund(x=solbdat,y=treatment,type="bx",
             main="BoxPlot of The 12 Highest Abundant Species", 
             xlab="Abundance [sqrt scale]",ylab="",
             transformation="sqrt",t.lab="o",shift=TRUE)

}
\keyword{hplot}

