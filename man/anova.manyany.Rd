\name{anova.manyany}
\alias{anova.manyany}
\alias{print.anova.manyany}
\title{Analysis of Deviance for Many Univariate Models Fitted to Multivariate Abundance Data}
\description{Compute an analysis of deviance table for many univariate model fits. Slowly!}
\usage{
\method{anova}{manyany}(object, \dots, nBoot=99, p.uni="none", block=object1$block, nCores=1,
       bootID=NULL, replace=TRUE) 

\method{print}{anova.manyany}(x, ...) 
}
\arguments{
\item{object}{of class \code{manyany} under the null hypothesis, typically the result of a call to \code{\link{manyany}}.}
\item{\dots}{other generic \code{anova} methods. NEEDS TO INCLUDE A SECOND \code{manyany} object for the alternative hypothesis to be tested.}
\item{nBoot}{the number of Bootstrap iterations, default is \code{nBoot=99}.}
\item{p.uni}{whether to calculate univariate test statistics and their P-values.\cr
"none" = No univariate P-values (default) \cr
"unadjusted" = A test statistic and (ordinary unadjusted) P-value are reported
for each response variable. If the \code{manyany} object is compositional (\code{composition=TRUE}), this option is unavailable as yet.
}
\item{block}{a factor specifying the sampling level to be resampled. Default is resampling rows (if
composition=TRUE in the manyany command, this means resampling rows of data as originally sent to manyany).}
\item{nCores}{Number of cores to use for computations (for parallel computing).}
\item{bootID}{A user-entered matrix of indices for which observations to use in which resample. Bootstrap
resamples in rows, observations in columns. When specified, overwrites \code{nBoot} and \code{block}. Default is NULL.}
\item{replace}{whether to sample with or without replacement, as in the \code{\link{sample}} function.
\code{= FALSE} for PIT-permutation, \code{ = TRUE} for PIT-trap.}
\item{x}{\code{anova.manyany} object to be printed.}
}

\details{
The \code{anova.manyany} function returns a table summarising the statistical significance
of a fitted manyany model under the alternative hypothesis (\code{object2}) as compared
to a fit under the null hypothesis (\code{object}). Typically the alternative model
is nested in the null although it doesn't need to be (but consider seriously if what
you are doing makes sense if they are not nested).

This function is quite computationally intensive, and a little fussy - it is an early
version we hope to improve on. Feedback welcome!

This function behaves a lot like \code{\link{anova.manyglm}}, the most conspicuous differences
being in flexibility and computation time. Since this function is based on \code{manyany},
it offers much greater flexibility in terms of types of models that can be fitted (most
fixed effects model with \code{predict} and \code{family} arguments could be accommodated).
For information on the different types of data that can be modelled using manyany, see
\code{\link{manyany}}.

However this flexibility comes at considerable cost in terms of computation time, and the
default \code{nBoot} has been set to 99 to reflect this (although rerunning at 999 is
recommended). Other more cosmetic differences from \code{\link{anova.manyglm}} are that
two and only two models can be supplied as input here; adjusted univariate P-values
are not yet implemented; and the range of test statistics and resampling algorithms is
more limited. All test statistics constructed here are sum-of-likelihood ratio statistics
as in Warton et al (2012), and the resampling method used here is the PIT-trap (short
for 'probability integral transform residual bootstrap', Warton et al 2017).

To check model assumptions, use \code{\link{plot.manyany}}.

The \code{block} argument allows for block resampling, such that valid inferences can
be made across independent blocks of correlated sets of observations.
For example, if data have multiple rows of records for each site, e.g. multi-species
data with entries for different species on different rows, you can use your site ID
variable as the block argument to resample sites, for valid cross-site inferences despite
within-site species correlation. Well, valid assuming sites are independent. You could
do similarly for a repeated measures design to make inferences robust to temporal autocorrelation.
Note that \code{block} needs to be balanced, e.g. equal number of species entries for
each site (i.e. include rows for zero abundances too). 

The \code{anova.manyany} function is designed specifically for high-dimensional data
(that, is when the number of variables p is not small compared to the number of observations
N). In such instances a correlation matrix is computationally intensive to estimate
and is numerically unstable, so by default the test statistic is calculated assuming
independence of variables. Note however that the resampling scheme used ensures that
the P-values are approximately correct even when the independence assumption is not
satisfied. 

Rather than stopping after testing for multivariate effects, it is often of interest
to find out which response variables express significant effects. Univariate statistics
are required to answer this question, and these are reported if requested. Setting \code{p.uni="unadjusted"}
returns resampling-based univariate P-values for all effects as well as the multivariate
P-values, if \code{composition=FALSE}. There are currently no univariate P-value options
when \code{composition=TRUE} (it's not entirely clear how such P-values should be obtained)
and if univariate P's are of interest why not rerun the model with \code{composition=FALSE}.

A current limitation of the function is that \code{composition} needs to be set to
the same value in each manyany object being compared - it is not currently possible
to compare models with and without a compositional term in them.
}

\value{
\item{stat}{the observed value of the test statistic.}
\item{p}{the P-value as estimated from \code{nBoot} resamples.}
\item{stat.i}{the values of the test statistic in each of the \code{nBoot} resamples.}
\item{p.i}{the P-value in each of the \code{nBoot} resamples.}
\item{p.uni}{the \code{p.uni} argument supplied.}

If \code{p.uni="unadjusted"} the output list also contains\cr
\item{uni.test}{a table showing the test statistics of the univariate tests.}
\item{uni.p}{a table showing the p-values of the univariate tests.}
\item{statj.i}{a matrix of values of the univariate test statistics in each of the \code{nBoot} resamples.}
}

\section{Warning}{
The comparison between two or more models by \code{anova.manyglm} will only be valid if they are fitted to the same dataset. This may be a problem if there are missing values and R's default of \code{na.action = na.omit} is used. 
}

\references{

Warton D. I., Wright S., and Wang, Y. (2012). Distance-based multivariate analyses confound location and dispersion effects. \emph{Methods in Ecology and Evolution}, 3(1), 89-101.

Warton D. I., Thibaut L., Wang Y. A.  (2017). The PIT-trap - A "model-free" bootstrap procedure for inference about regression models with discrete, multivariate responses. \emph{PLoS One}, 12(7), e0181790.

}
\author{
David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{manyany}}, \code{\link[mvabund:anova.manyglm]{anova.manyglm}}.
}
\examples{
## Try fitting Tikus Islands data with Tweedie models with power parameter 1.5,
## to test for compositional effect:
data(tikus)
coral <- as.matrix(tikus$abund[1:20,])
sumSpp = apply(coral>0,2,sum)

coral <- coral[,sumSpp>6] ## cutting to just species with seven(!) or more presences to cut
## computation time.  Maybe rerun with less (e.g. 4 or more presences) if curious and patient.
coralX <- tikus$x[1:20,]

require(tweedie)
require(statmod)

ftTimeRep <- manyany(coral ~ time+rep, "glm", data=coralX, 
  family=tweedie(var.power=1.5, link.power=0), var.power=1.5, composition=TRUE)

ftRep <- manyany(coral ~ rep, "glm", data=coralX, 
  family=tweedie(var.power=1.5, link.power=0), var.power=1.5, composition=TRUE)
anova(ftRep,ftTimeRep,nBoot=9) #this takes a few seconds to run even for just 9 resamples
## This should be rerun for nBoot=999, which could a few minutes...

\dontrun{library(ordinal)}
## First construct an ordinal dataset:
\dontrun{spidOrd = spider$abund}
\dontrun{spidOrd[spider$abund>1 & spider$abund<=10]=2}
\dontrun{spidOrd[spider$abund>10]=3}
## Now fit a model using the clm function:
\dontrun{manyOrd=manyany(spidOrd~bare.sand+fallen.leaves,"clm",data=spider$x)}
## Test to see if fallen.leaves needs to be there:
\dontrun{manyOrd0=manyany(spidOrd~bare.sand,"clm",data=spider$x)}
\dontrun{anova(manyOrd0,manyOrd,nBoot=19)}

}
\keyword{models}
\keyword{regression}
\keyword{multivariate}
