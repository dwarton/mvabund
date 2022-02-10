\name{summary.manylm}
\alias{summary.manylm}
\alias{print.summary.manylm}
\title{Summarizing Linear Model Fits for Multivariate Abundance Data}
\description{\code{summary} method for class "manylm" - computes a table 
summarising the statistical significance of different multivariate terms 
in a linear model fitted to high-dimensional data, such as multivariate 
abundance data in ecology.}

\usage{
\method{summary}{manylm}(object, nBoot=999, resamp="residual", 
   test="F", cor.type=object$cor.type, block=NULL, shrink.param=NULL, 
   p.uni="none", studentize=TRUE, R2="h", show.cor = FALSE, 
   show.est=FALSE, show.residuals=FALSE, symbolic.cor=FALSE, 
   rep.seed=FALSE, tol=1.0e-6, \dots )
  
\method{print}{summary.manylm}(
   x, digits = max(getOption("digits") - 3, 3), 
   signif.stars=getOption("show.signif.stars"), 
   dig.tst=max(1, min(5, digits - 1)), 
   eps.Pvalue=.Machine$double.eps, \dots )
}

\arguments{
\item{object}{an object of class "manylm", usually, a result of a call to
\code{\link{manylm}}.}
\item{nBoot}{the number of Bootstrap iterations, default is \code{nBoot=999}.}
\item{resamp}{the method of resampling used. Can be one of "case" (not yet available),"residual" (default), "perm.resid", "score" or "none". See Details.}
\item{test}{ the test to be used. Possible values are: "LR" = likelihood ratio statistic (default) and "F" = Lawley-Hotelling trace statistic.\cr
Note that if all variables are assumed independent (\code{cor.shrink="I"}) then "LR" is equivalent to LR-IND and "F" is the sum-of-F statistics from Warton & Hudson (2004).}
\item{cor.type}{ structure imposed on the estimated correlation matrix under the fitted model. Can be "I"(default), "shrink", or "R". See Details.}
\item{block}{A factor specifying the sampling level to be resampled. Default is resampling rows.}
\item{shrink.param}{shrinkage parameter to be used if \code{cor.type="shrink"}. If not supplied, but needed, it will be estimated from the data by Cross
Validation using the normal likelihood as in Warton (2008).}
\item{p.uni}{whether to calculate univariate test statistics and their P-values, and if so, what type. \cr
"none" = no univariate P-values (default) \cr
"unadjusted" = a test statistic and (ordinary unadjusted) P-value is reported
for each response variable. \cr
"adjusted" = Univariate P-values are adjusted for multiple testing, using a
step-down resampling procedure.
}
\item{studentize}{logical, whether studentized residuals or residuals should beused for simulation in the resampling steps. This option is not used in case resampling. }
\item{R2}{
 the type of R^2 (correlation coefficient) that should be shown, can be one of: \cr
 "h" = Hooper's R^2 = tr(SST^(-1)SSR)/p   \cr
 "v" = vector R^2 = det(SSR)/det(SST)\cr
 "n" = none
}
\item{show.cor}{
logical, if \code{TRUE}, the correlation matrix of the estimated parameters is
returned and printed.}
\item{show.est}{logical. Whether to show the estimated model parameters.}
\item{show.residuals}{logical. Whether to show residuals/a residual summary.}
\item{symbolic.cor}{logical. If \code{TRUE}, print the correlations in a symbolic form rather than as numbers.}
\item{rep.seed}{logical. Whether to fix random seed in resampling data. Useful for simulation or diagnostic purposes.}
\item{tol}{the tolerance used in estimations.}

\item{x}{an object of class \code{"summary.manylm"}, usually, a result of a call to \code{summary.manylm}.}
\item{digits}{the number of significant digits to use when printing.}
\item{signif.stars}{logical. If \code{TRUE}, \sQuote{significance stars} are printed for each coefficient.}
\item{dig.tst}{the number of digits to round the estimates of the model parameters.}
\item{eps.Pvalue}{a numerical tolerance for the formatting of p values.}
\item{\dots}{for \code{summary.manyglm} method, these are additional arguments including: \cr
\code{bootID} - this matrix should be integer numbers where each row specifies bootstrap id's in each resampling run. When \code{bootID} is supplied, \code{nBoot} is set to the number of rows in \code{bootID}. Default is \code{NULL}.\cr
for \code{print.summary.manyglm} method, these are optional further arguments passed to or from other methods. See \code{\link{print.summary.glm}} for more details.}
}

\details{ 
The \code{summary.manylm} function returns a table summarising the statistical 
significance of each multivariate term specified in the fitted manylm model. 
For each model term, it returns a test statistic as determined by the argument 
\code{test}, and a P-value calculated by resampling rows of the data using a 
method determined by the argument \code{resamp}. The four possible resampling methods are residual-permutation (Anderson and Robinson (2001)), score resampling (Wu (1986)), case and residual resampling (Davison and Hinkley (1997, chapter 6)), and involve resampling under the alternative hypothesis. These methods ensure approximately valid inference even when the correlation between variables has been misspecified, and for case and score resampling, even when the equal variance assumption of linear models is invalid. By default, studentized residuals (r_i/sqrt(1-h_ii)) are used in residual and score resampling, although raw residuals could be used via the argument \code{studentize=FALSE}. If \code{resamp="none"}, p-values cannot be calculated, however the test statistics are returned.

If you have a specific hypothesis of primary interest that you want to test,
then you should use the \code{anova.manylm} function, which can resample rows
of the data under the null hypothesis and so usually achieves a better
approximation to the true significance level.

To check model assumptions, use \code{plot.manylm}.

The \code{summary.manylm} function is designed specifically for high-dimensional data (that, is when the number of variables p is not small compared to the number of observations N). In such instances a correlation matrix is computationally intensive to estimate and is numerically unstable, so by default the test statistic is calculated assuming independence of variables (\code{cor.type="I"}). Note however that the resampling scheme used ensures that the P-values are approximately correct even when the independence assumption is not satisfied. However if it is computationally feasible for your dataset, it is recommended that you use \code{cor.type="shrink"} to account for correlation between variables, or \code{cor.type="R"} when p is small. The \code{cor.type="R"} option uses the unstructured correlation matrix (only possible when N>p), such that the standard classical multivariate test statistics are obtained. Note however that such statistics are typically numerically unstable and have low power when p is not small compared to N. The \code{cor.type="shrink"} option applies ridge regularisation (Warton 2008), shrinking the sample correlation matrix towards the identity, which improves its stability when p is not small compared to N. This provides a compromise between \code{"R"} and \code{"I"}, allowing us to account for correlation between variables, while using a numerically stable test statistic that has good properties. The shrinkage parameter by default is estimated by cross-validation using the multivariate normal likelihood function, although it can be specified via \code{shrink.param} as any value between 0 and 1 (0="I" and 1="R", values closer towards 0 indicate more shrinkage towards "I"). The validation groups are chosen by random assignment and so you may observe some slight variation in the estimated shrinkage parameter in repeat analyses. See \code{\link{ridgeParamEst}} for more details.

Rather than stopping after testing for multivariate effects, it is often of interest to find out which response variables express significant effects. Univariate statistics are required to answer this question, and these are reported if requested. Setting \code{p.uni="unadjusted"} returns resampling-based univariate P-values for all effects as well as the multivariate P-values, whereas  \code{p.uni="adjusted"} returns adjusted P-values (that have been adjusted for multiple testing), calculated using a step-down resampling algorithm as in Westfall & Young (1993, Algorithm 2.8). This method provides strong control of family-wise error rates, and makes use of resampling (using the method controlled by \code{resample}) to ensure inferences take into account correlation between variables.

A multivariate R^2 value is returned in output, but there are many ways to define a multivariate R^2. The type of R^2 used is controlled by the \code{R2} argument. If \code{cor.shrink="I"} then all variables are assumed independent, a special case in which Hooper's R^2 returns the average of all univariate R^2 values, whereas the vector R^2 returns their product. 

\code{print.summary.manylm} tries to be smart about formatting the coefficients, \code{genVar}, etc. and additionally gives \sQuote{significance stars} if
\code{signif.stars} is \code{TRUE}.
}

\value{
summary.manylm returns an object of class "summary.manyglm", a list with components 
\item{call}{the component from \code{object}.}
\item{terms}{the terms object used.}
\item{show.residuals}{the supplied argument.}
\item{show.est}{the supplied argument.}
\item{p.uni}{the supplied argument.}
\item{test}{the supplied argument.}
\item{cor.type}{the supplied argument.}
\item{resample}{the supplied argument.}
\item{nBoot}{the supplied argument.}
\item{rankX}{the rank of the design matrix}
\item{residuals}{the model residuals}
\item{genVar}{the estimated generalised variance}
\item{est}{the estimated model coefficients}
\item{shrink.param}{
the shrinkage parameter. Either the value of the argument with the same name or
if this was not supplied the estimated shrinkage parameter.
}
\item{aliased}{named logical vector showing if the original coefficients are
aliased.} 
\item{df}{ a 3-vector of the rank of the model and the number of residual degrees of freedom, plus number of non-aliased coefficients.}

If the argument \code{test} is not \code{NULL} then the list also
included the components
\item{coefficients}{a matrix containing the test statistics and the p-values.}
\item{n.iter.sing}{the number of iterations that were skipped due to singularity of the design matrix caused by case resampling.}

If furthermore the Design matrix is neither empty nor consists of the Intercept only, the following adddional components are included:
\item{r.squared}{the calculated correlation coefficient.}
\item{R2}{a character that describes which type of correlation coefficient was
calculated.}
\item{statistic}{a matrix containing the results of the overall test.}
\item{cov.unscaled}{
the unscaled (\code{dispersion = 1}) estimated covariance matrix of the
estimated coefficients.
}

If the argument \code{show.cor} is \code{TRUE} the following adddional
components are returned:
\item{correlation}{the (p*q) by (p*q) correlation matrix, with p being the
number of columns of the design matrix and q being the number of response
variables. Note that this matrix can be very big.}
}

\references{
Anderson, M.J. and J. Robinson (2001).
Permutation tests for linear models.
\emph{Australian and New Zealand Journal of Statistics} 43, 75-88.

Davison, A. C. and Hinkley, D. V. (1997)
\emph{Bootstrap Methods and their Application.}
Cambridge University Press, Cambridge.

Warton D.I. (2008). Penalized normal likelihood and ridge regularization of 
correlation and covariance matrices. \emph{Journal of the American Statistical 
Association} 103, 340-349.

Warton D.I. and Hudson H.M. (2004). A MANOVA statistic is just as powerful as
distance-based statistics, for multivariate abundances. 
\emph{Ecology} 85(3), 858-874.

Westfall, P. H. and Young, S. S. (1993) \emph{Resampling-based multiple 
testing.} John Wiley & Sons, New York.

Wu, C. F. J. (1986) Jackknife, Bootstrap and Other Resampling Methods in
Regression Analysis.
\emph{The Annals of Statistics} 14:4, 1261-1295.
}
\author{
Yi Wang, Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{manylm}}, \code{\link{anova.manylm}}, \code{\link{plot.manylm}}
}
\examples{
data(spider)
spiddat <- log(spider$abund+1)

## Estimate the coefficients of a multivariate linear model:
fit <- manylm(spiddat~., data=spider$x)

## To summarise this multivariate fit, using score resampling to
## and F Test statistic to estimate significance:
summary(fit, resamp="score", test="F")

## Instead using residual permutation with 2000 iteration, using the sum of F 
## statistics to estimate multivariate significance, but also reporting 
## univariate statistics with adjusted P-values:
summary(fit, resamp="perm.resid", nBoot=2000, test="F", p.uni="adjusted")

## Obtain a summary of test statistics using residual resampling, accounting 
## for correlation between variables but shrinking the correlation matrix to 
## improve its stability and showing univariate p-values:
summary(fit, cor.type="shrink", p.uni="adjusted")

}
\keyword{models}
\keyword{regression}
\keyword{multivariate}

