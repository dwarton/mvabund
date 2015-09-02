\name{plot.manylm}
\alias{plot.manylm}
\alias{plot.manyglm}
\title{Plot Diagnostics for a manylm or a manyglm Object}
\description{
  Four plots (selectable by \code{which}) are currently available: a plot
  of residuals against fitted values, a Normal Q-Q plot, 
  a Scale-Location plot of \eqn{\sqrt{| residuals |}} against fitted values, 
  a plot of Cook's distances versus row labels.  
  By default, all of them are provided. 

  The function is not yet available for manyglm object
}
\usage{
\method{plot}{manylm}(
  x, res.type="pearson", which=1:4, caption=c("Residuals vs Fitted", "Normal Q-Q", 
  "Scale-Location", "Cook's distance"), overlay=TRUE, 
  n.vars=Inf, var.subset=NULL, sub.caption=NULL, studentized= TRUE, ...)

\method{plot}{manyglm}(
  x, res.type="pit.norm", which=1, caption=c("Residuals vs Fitted", "Normal Q-Q", 
  "Scale-Location", "Cook's distance"), overlay=TRUE, 
  n.vars=Inf, var.subset=NULL, sub.caption=NULL, ...)

}
\arguments{
 \item{x}{ \code{manylm} object or \code{manyglm} object, typically the 
  result of a call to \code{\link{manylm}} or \code{manyglm}.}
 \item{res.type}{type of residuals to plot. By default, \code{res.type="pit-norm"} uses Dunn-Smyth
 residuals (Dunn & Smyth 2996), related to the probability integral transform, for manyglm.
 These residuals are especially recommended for presence-absence data or discrete data. }
  \item{which}{if a subset of the plots is required, specify a subset of
    the numbers \code{1:4}.}
  \item{caption}{captions to appear above the plots}
  \item{overlay}{logical, whether or not the different variables should be
    overlaid on a single plot.}
  \item{n.vars}{the number of variables to include in the plot.}
  \item{var.subset}{the variables to include in the plot.}
  \item{sub.caption}{common title---above figures if there are multiple;
    used as \code{sub} (s.\code{title}) otherwise.  If \code{NULL},
    as by default, a possible shortened version of
    \code{deparse(x$call)} is used.}
  \item{\dots}{other parameters to be passed through to plotting
    functions.}
  \item{studentized}{logical indicating whether studentized or standardized residuals should be used for plot 2 and 3.}

}
\details{
  \code{plot.manylm} is used to check the linear model assumptions that are made
  when fitting a model via \code{manylm}. Similarly, \code{plot.manyglm} checks
  the generalised linear model assumptions made when using \code{manyglm}.
  As in Wang et al (2012), you should check the residual vs fits plot for no pattern
  (hence no suggestion of failure of key linearity and mean-variance assumptions).
  For manylm fits of small datasets, it is desirable that residuals on the normal Q-Q plot be close
  to a straight line, although in practice the most important thing is to make
  sure there are no big outliers and no suggestion of strong skew in the data.
  
  The recommended \code{res.type} option for manyglm calls, "pit-norm", uses randomised quantile or "Dunn-Smyth"
  residuals (Dunn & Smyth 1996). Note that for discrete data, these residuals
  involve random number generation, and will not return identical results on replicate runs - so it is recommended
  that you plot your data a few times to check if any pattern shows up consistently across replicate plots.
  The other main residual option is "pearson", Pearson residuals. Note that all res.type options
  are equivalent for manylm.

  Some technical details on usage of this function:  \cr
   \code{sub.caption} - by default the function call - is shown as
  a subtitle (under the x-axis title) on each plot when plots are on
  separate pages, or as a subtitle in the outer margin (if any) when
  there are multiple plots per page.  % check that
  
  The \sQuote{Scale-Location} plot, also called \sQuote{Spread-Location} or
  \sQuote{S-L} plot, takes the square root of the absolute residuals in
  order to diminish skewness (\eqn{\sqrt{| E |}} is much less skewed
  than \eqn{| E |} for Gaussian zero-mean \eqn{E}).
  
  If \code{studentized=FALSE} the \sQuote{S-L}, the Q-Q, and the Residual-Leverage 
  plot, use \emph{standardized} residuals which have identical variance 
  (under the hypothesis) otherwise \emph{studentized} residuals are used.  

 Unlike other plotting functions \code{plot.manylm} and \code{plot.manyglm} 
 respectively do not have a subset argument, the subset needs to be specified 
 in the \code{manylm} or respectively \code{manyglm} function.

 For all arguments that are formally located after the position of \code{\dots},
 positional matching does not work.

 For restrictions on \code{filename} see R's help on eps/pdf/jpeg.
 Note that \code{keep.window} will be ignored if \code{write.plot} is 
 not \code{show}.
}
\references{

Dunn, P.K., & Smyth, G.K. (1996). Randomized quantile residuals. Journal of Computational and
Graphical Statistics 5, 236-244.

Wang Y., Naumann U., Wright S.T. & Warton D.I. (2012). mvabund - an R package for model-based 
analysis of multivariate abundance data. Methods in Ecology and Evolution 3, 471-474.

}
\author{
Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{manylm}}
}
\examples{
require(graphics)

data(spider)
spiddat <- mvabund(spider$abund)
X <- spider$x

## plot the diagnostics for the linear fit of the spider data
spidlm <- manylm(spiddat~X)
plot(spidlm,which=1:2,col.main="red",cex=3,overlay=FALSE)

plot(spidlm,which=1:4,col.main="red",cex=3,overlay=TRUE)

## plot the diagnostics for Poisson and negative binomial regression of the spider data
glmP.spid <- manyglm(spiddat~X, family="poisson")
plot(glmP.spid, which=1) #note the marked fan-shape on the plot

glmNB.spid <- manyglm(spiddat~X, family="negative.binomial")
plot(glmNB.spid, which=1) #no fan-shape
plot(glmNB.spid, which=1) #note the residuals change on re-plotting, but no consistent trend

}
\keyword{hplot}
\keyword{regression}
\keyword{multivariate}
