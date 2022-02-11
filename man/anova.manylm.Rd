\name{anova.manylm}
\alias{anova.manylm}
\alias{print.anova.manylm}
\title{ANOVA for Linear Model Fits for Multivariate Abundance Data}
\description{
\code{anova} method for class "manylm" - computes an analysis of variance 
table for one or more linear model fits to high-dimensional data, such as 
multivariate abundance data in ecology.
}

\usage{
\method{anova}{manylm}(object, \dots, resamp="perm.resid", test="F", p.uni="none",
        nBoot=999, cor.type=object$cor.type,
        block=NULL, shrink.param=object$shrink.param, 
	studentize=TRUE, calc.rss = FALSE, tol=1.0e-10, rep.seed=FALSE, bootID=NULL )
\method{print}{anova.manylm}(
    x, digits = max(getOption("digits") - 3, 3),
    signif.stars = getOption("show.signif.stars"),
    dig.tst = max(1, min(5, digits - 1)),
    eps.Pvalue = .Machine$double.eps, \dots)
}
\arguments{
\item{object}{
object of class \code{manylm}, usually, a result of a call to \code{\link{manylm}}.
}
\item{\dots}{
for the \code{anova.manylm} method, these are optional further objects of class \code{manylm}, which are usually a result of a call to \code{\link{manylm}}.
for the \code{print.anova.manylm} method these are optional further arguments
passed to or from other methods.
}
\item{nBoot}{the number of iterations in resampling. Default is 999 for P-values as fractions out of 1000.}
\item{resamp}{the method of resampling used. Can be one of "perm.resid" (default), "residual", "score", "case". See Details.}
\item{test}{the test to be used. Possible values are: \code{"LR"} = likelihood ratio statistic, \code{"F"} = Lawley-Hotelling trace statistic or \code{NULL} for no test.}
\item{cor.type}{structure imposed on the estimated correlation matrix under the fitted model. Can be "I"(default), "shrink", or "R". See Details.}
\item{block}{A factor specifying the sampling level to be resampled. Default is resampling rows.}
\item{shrink.param}{shrinkage parameter to be used if \code{cor.type="shrink"}. If not supplied, but needed, it will be estimated by estimated from the data by Cross Validation using the normal likelihood as in Warton (2008).}
\item{p.uni}{whether to calculate univariate test statistics and their P-values, and if so, what type. \cr
"none" = no univariate P-values (default) \cr
"unadjusted" = A test statistic and (ordinary unadjusted) P-value is reported
for each response variable. \cr
"adjusted" = Univariate P-values are adjusted for multiple testing, using a step-down resampling procedure.
}
\item{studentize}{logical. Whether studentized residuals should be used to simulate the data in the resampling steps. This option is not used in case resampling.}
\item{calc.rss}{logical. Whether the Residual Sum of Squares should be calculated.}
\item{tol}{the sensitivity in calculations near 0.}
\item{rep.seed}{logical. Whether to fix random seed in resampling data. Useful for simulation or diagnostic purposes.}
\item{bootID}{an integer matrix where each row specifies bootstrap id's in each resampling run. When \code{bootID} is supplied, \code{nBoot} is set to the number of rows in \code{bootID}. Default is \code{NULL}.}
\item{x}{an object of class \code{"anova.manylm"}, usually, a result of a call to \code{anova.manylm}.}
\item{digits}{the number of significant digits to use when printing.}
\item{signif.stars}{logical. If \code{TRUE}, \sQuote{significance stars} are
printed for each coefficient.}
\item{dig.tst}{the number of digits to round the estimates of the model parameters.}
\item{eps.Pvalue}{a numerical tolerance for the formatting of p values.}
}

\details{
The \code{anova.manylm} function returns a table summarising the statistical significance of a fitted manylm model, or of the differences between several nested models fitted to the same dataset. If one model is specified, sequential test statistics (and P values) are returned for that fit. If more than one object is specified, the table contains test statistics (and P values) comparing their fits.

The test statistics are determined by the argument \code{test}, and the P-values are calculated by resampling rows of the data using a method determined by the argument \code{resamp}. The four possible resampling methods are residual-permutation (Anderson and Robinson (2001)), score resampling (Wu (1986)), case and residual resampling (Davison and Hinkley (1997, chapter 6)), and involve resampling under the null hypothesis (except for case resampling). These methods ensure approximately valid inference even when the correlation between variables has been misspecified, and for case and score resampling, even when the equal variance assumption of linear models is invalid. By default, studentised residuals (r_i/sqrt(1-h_ii)) are used in residual and score resampling, although raw residuals could be used via the argument \code{studentize=FALSE}. 
If \code{resamp="none"}, p-values cannot be calculated, however the test statistics are returned.

If you do not have a specific hypothesis of primary interest that you want to test, and are instead interested in which model terms are statistically significant, then the \code{summary.manylm} function is more appropriate.

More than one object should only be specified when the models are nested. In this case the ANOVA table has a column for the residual degrees of freedom and a column for change in degrees of freedom. It is conventional to list the models from smallest to largest, but this is up to the user.

To check model assumptions, use \code{plot.manylm}.

The \code{anova.manylm} function is designed specifically for high-dimensional data (that, is when the number of variables p is not small compared to the number of observations N). In such instances a correlation matrix is computationally intensive to estimate and is numerically unstable, so by default the test statistic is calculated assuming independence of variables (\code{cor.type="I"}). Note however that the resampling scheme used ensures that the P-values are approximately correct even when the independence assumption is not satisfied. However if it is computationally feasible for your dataset, it is recommended that you use \code{cor.type="shrink"} to account for correlation between variables, or \code{cor.type="R"} when p is small. The \code{cor.type="R"} option uses the unstructured correlation matrix (only possible when N>p), such that the standard classical multivariate test statistics are obtained. Note however that such statistics are typically numerically unstable and have low power when p is not small compared to N. The \code{cor.type="shrink"} option applies ridge regularisation (Warton 2008), shrinking the sample correlation matrix towards the identity, which improves its stability when p is not small compared to N. This provides a compromise between \code{"R"} and \code{"I"}, allowing us to account for 
correlation between variables, while using a numerically stable test statistic that has good properties. The shrinkage parameter by default is estimated by cross-validation using the multivariate normal likelihood function, although it can be specified via \code{shrink.param} as any value between 0 and 1 (0="I" and 1="R", values closer towards 0 indicate more shrinkage towards "I"). The validation groups are chosen by random assignment and so you may observe some slight variation in the estimated shrinkage parameter in repeat analyses. 
See \code{\link{ridgeParamEst}} for more details.

Rather than stopping after testing for multivariate effects, it is often of interest to find out which response variables express significant effects. Univariate statistics are required to answer this question, and these are reported if requested. Setting \code{p.uni="unadjusted"} returns the resampling-based univariate ANOVA P-values as well as the multivariate P-values, whereas  \code{p.uni="adjusted"} returns adjusted ANOVA P-values (that have been adjusted for 
multiple testing), calculated using a step-down resampling algorithm as in Westfall & Young (1993, Algorithm 2.8). This method provides strong control of family-wise error rates, and makes use of resampling (using the method controlled by \code{resampling}) to ensure inferences take into account correlation between variables.

}

\value{
An object of class \code{"anova.manylm"}. A list containing at least:  \cr
\item{p.uni}{the supplied argument.}
\item{test}{the supplied argument.}
\item{cor.type}{the supplied argument.}
\item{resample}{the supplied argument.}
\item{nBoot}{the supplied argument.}
\item{calc.rss}{the supplied argument.}
\item{table}{the data frame containing the anova table.}
\item{shrink.param}{the supplied argument.}
\item{n.bootsdone}{the number of bootstrapping iterations that were done,
i.e. had no error.}
\item{n.iter.sing}{the number of bootstrap iterations where the resampled
design matrix was singular and could only be used partly.}
\item{one}{
logical. whether the anova table was calculated for one \code{manylm} object
or for several \code{manylm} objects.
}

If \code{p.uni="adjusted"} or \code{p.uni="unadjusted"} the output list also
contains: \cr
\item{uni.test}{a table showing the test statistics of the univariate tests}
\item{uni.p}{a table showing the p-values of the univariate tests}

The print method for \code{anova.manylm} objects prints the output in a
\sQuote{pretty} form.
}

\references{
Anderson, M.J. and J. Robinson (2001).
Permutation tests for linear models.
\emph{Australian and New Zealand Journal of Statistics} 43, 75-88.

Davison, A. C. and Hinkley, D. V. (1997)
\emph{Bootstrap Methods and their Application.}
Cambridge University Press, Cambridge.

Warton D.I. (2008). Penalized normal likelihood and ridge regularization of 
correlation and covariance matrices. \emph{Journal of the American 
Statistical Association} 103, 340-349.

Warton D.I. and Hudson H.M. (2004). A MANOVA statistic is just as powerful as
distance-based statistics, for multivariate abundances. 
\emph{Ecology} 85(3), 858-874.

Westfall, P. H. and Young, S. S. (1993) 
\emph{Resampling-based multiple testing.} John Wiley & Sons, New York.

Wu, C. F. J. (1986) Jackknife, Bootstrap and Other Resampling Methods in
Regression Analysis.
\emph{The Annals of Statistics} 14:4, 1261-1295.
}
\author{
Yi Wang, Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{manylm}}, \code{\link{summary.manylm}}, \code{\link{plot.manylm}}
}

\examples{
## Load the spider dataset:
data(spider)
spiddat <- log(spider$abund+1)
spiddat <- mvabund(spiddat)
spidx <- as.matrix(spider$x)

## Fit several multivariate linear models:
fit <- manylm( spiddat ~ spidx ) # model with all explanatory variables

## Use the default residual resampling to test the significance of this model:
## return summary of the manylm model
anova(fit)

## intercept model
fit0 <- manylm(spiddat ~ 1)
## include soil and leaf variables
fit1 <- update(fit0, . ~ . + spidx[, c(1, 3)])
## include moss variables
fit2 <- update(fit1, . ~ . + spidx[, 4]) 

## Use (residual) resampling to test the significance of these models, 
## accounting for correlation between variables by shrinking 
## the correlation matrix to improve its stability:
anova(fit, fit0, fit1, fit2, cor.type="shrink")

## Use the sum of F statistics to estimate multivariate significance from 
## 4999 resamples, and also reporting univariate statistics with 
## adjusted P-values:
anova(fit, fit0, fit1, fit2, nBoot=4999, test="F", p.uni="adjust")

}
\keyword{regression}
\keyword{models}
\keyword{multivariate}

