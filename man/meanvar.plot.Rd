\name{meanvar.plot}
\alias{meanvar.plot}
\alias{meanvar.plot.mvabund}
\alias{meanvar.plot.mvformula}

\alias{meanvar.plot,mvabund-method}
\alias{meanvar.plot,mvformula-method}
\alias{meanvar.plot,formula-method}
\alias{meanvar.plot,matrix-method}
\alias{meanvar.plot,data.frame-method}

\title{
Construct Mean-Variance plots for Multivariate Abundance Data
}
\description{
Construct mean-variance plots, separately for each column of the input
 data, and separately for each level of any input factor that is 
 given (via a formula). This function was specially written for high 
 dimensional data where there are many correlated variables exhibiting 
 a mean-variance structure, in particular, multivariate abundance 
 data in ecology.
}
\usage{
meanvar.plot(x, \dots)

\method{meanvar.plot}{mvabund}(
  x, n.vars=NULL, var.subset=NULL, subset=NULL, table=FALSE, ...)

\method{meanvar.plot}{mvformula}(
  x, n.vars = NULL, var.subset=NULL, subset=NULL, table=FALSE, 
  overall.main=NULL, overlay=TRUE, ...) 
}

\arguments{
\item{x}{an mvabund objects or a Model Formula (can be a formula or a
mvformula) to be used.}
\item{n.vars}{the number of variables to include in the plot.}
\item{var.subset}{vector of indices indicating the variables to be included on the plot, (default: the \code{n.vars} most abundant variables).}
\item{subset}{an optional vector specifying a subset of observations to be used.}
\item{table}{logical, whether a table of the Means and Variances should be returned}
\item{overall.main}{an overall title for the window.}
\item{overlay}{logical, whether overall means/variances for all variables are calculated and drawn on a single plot or calculated and plotted separately for different variables.}
\item{\dots}{arguments to be passed to or from other methods.} 
}
\details{
\code{meanvar.plot} calculates a mean-variance plot for a dataset with many variables (e.g., Warton D. I., Wright S., and Wang, Y. (2012)).

The mean values and variances are calculated across all observations, unless a 
formula is given as the first argument which specifies a factor as the dependent 
variable. In this latter case the means and variances are calculated separately within the groups defined by these factors.

By default the means and variances of all variables (and all factor levels) are displayed on the same plot. If a formula is given and the explanatory variables contain factor variables, the mean values and variances for each factor level can be calculated and displayed either for all variables together or for each variable separately.

For the latter, set \code{overlay} to \code{FALSE}. The mean-variances corresponding to the different factors will be drawn in different colors, that can be chosen by specifying \code{col}. \code{col} can then either be a single color value (see \code{par}) with the number of values being at least the maximum number of levels of the factors. The same applies to \code{pch}.

If \code{mfrow} is \code{NULL} and \code{mfcol} is \code{NULL}, par("mfrow") is used. If \code{all.labels = FALSE}, only the x-axis labels at the bottom plot and the y-axis labels of plots on the right side of the window are printed if furthermore \code{main=NULL} only the graphics on the top contain the full title, the other ones an abreviated one.

Note, that if a log-transformation is used for displaying the data, a specific mean-variance relation will not be visible in the plot, if either the calculated mean is zero and \code{log!="x"} or \code{log!="xy"} or if the calculated variance
is zero and \code{log!="y"} or \code{log!="xy"}.

By default the y/x ratio of the axis, specified by \code{asp}, will be set to \code{1} if \code{log!="xy"}. If the mean-variance relation is not displayed on a log scale and \code{overlay} is \code{FALSE}, it is most often not advisable to specify \code{asp}, as there might not be one choice of \code{asp} that is sensible for each of the plots.

}

\value{
If \code{table} is \code{TRUE} a table of the Means and Variances is returned.
Otherwise, only the plot(s) is/are drawn.
}

\references{

Warton D. I., Wright S., and Wang, Y. (2012). Distance-based multivariate analyses confound location and dispersion effects. \emph{Methods in Ecology and Evolution}, 3(1), 89-101.


Warton D.I. (2008). Raw data graphing: an informative but under-utilized tool for the analysis of multivariate abundances. \emph{Austral Ecology} 33(3), 290-300.
}

\author{
Ulrike Naumann, Stephen Wright and David Warton <David.Warton@unsw.edu.au>.
}

\seealso{
\code{\link{plot.mvabund}}
\code{\link{plot.mvformula}}.
}

\examples{
require(graphics)

## Load the tikus dataset:
data(tikus)
tikusdat <- mvabund(tikus$abund)
year <- tikus$x[,1]

## Plot mean-variance plot for a mvabund object with a log scale (default):
meanvar.plot(tikusdat) 	

## Again but without log-transformation of axes:
meanvar.plot(tikusdat,log="") 	

## A mean-variance plot, data organised by year, 
## for 1981 and 1983 only, as in Figure 7a of Warton (2008):
is81or83 <- year==81 | year==83
meanvar.plot(tikusdat~year, subset=is81or83, col=c(1,10))
}
\keyword{hplot}

