\name{summary.manyglm}
\alias{summary.manyglm}
\alias{print.summary.manyglm}
\title{Summarizing Multivariate Generalized Linear Model Fits for Abundance Data}
\description{\code{summary} method for class "manyglm". }
\usage{
\method{summary}{manyglm}(object, resamp="pit.trap", test="wald", 
        p.uni="none", nBoot=999, cor.type=object$cor.type, block=NULL,
        show.cor = FALSE, show.est=FALSE, show.residuals=FALSE,
        symbolic.cor = FALSE,
        rep.seed = FALSE, 
        show.time=FALSE, show.warning=FALSE,\dots) 
\method{print}{summary.manyglm}(x, \dots) } \arguments{ \item{object}{an 
object of class \code{manyglm}, typically the result of a call to 
\code{\link{manyglm}}. } \item{resamp}{the method of resampling used. Can be 
one of "case", "perm.resid", "montecarlo" or "pit.trap" (default). See 
Details.}
\item{test}{the test to be used. If \code{cor.type="I"}, this can be 
one of "wald" for a Wald-Test or "score" for a Score-Test or "LR" for a 
Likelihood-Ratio-Test, otherwise only "wald" and "score" is allowed. The 
default value is "LR". }
\item{p.uni}{whether to calculate univariate test 
statistics and their P-values, and if so, what type. This can be one of the 
following options. \cr "none" = No univariate P-values (default) \cr 
"unadjusted" = A test statistic and (ordinary unadjusted) P-value are 
reported for each response variable. \cr "adjusted" = Univariate P-values are 
adjusted for multiple testing, using a step-down resampling procedure. } 
\item{nBoot}{the number of Bootstrap iterations, default is 
\code{nBoot=999}.}
\item{cor.type}{structure imposed on the estimated correlation matrix under the fitted model. Can be "I"(default), "shrink", or "R". See Details.}
\item{block}{A factor specifying the sampling level to be resampled. Default is resampling rows.}
\item{show.cor, show.est, show.residuals}{logical, if \code{TRUE}, 
the correlation matrix of the estimated parameters, or the estimated model 
parameters, or the residual summary is shown.}

\item{symbolic.cor}{logical. If \code{TRUE}, the correlation is printed in a symbolic form (see 
\code{\link{symnum}} ) rather than in numerical format. } 
\item{rep.seed}{logical. Whether to fix random seed in resampling data. Useful for simulation or diagnostic purposes.}

\item{show.time}{Whether to display timing information for the resampling 
procedure: "none" shows none, "all" shows all timing information and "total" 
shows only the overall time taken for the tests. } 
\item{show.warning}{logical. Whether to display warnings in the operation 
procedure.}
\item{\dots}{for \code{summary.manyglm} method, these are 
additional arguments including: \cr \code{rep.seed} - logical. Whether to fix 
random seed in resampling data. Useful for simulation or diagnostic 
purposes.\cr \code{bootID} - this matrix should be integer numbers where each 
row specifies bootstrap id's in each resampling run. When \code{bootID} is 
supplied, \code{nBoot} is set to the number of rows in \code{bootID}. Default 
is \code{NULL}.\cr for \code{print.summary.manyglm} method, these are 
optional further arguments passed to or from other methods. See 
\code{\link{print.summary.glm}} for more details.} \item{x}{an object of 
class "summary.manyglm", usually, a result of a call to 
\code{summary.manyglm}.} } 

\details{ The \code{summary.manyglm} function returns a table summarising the 
statistical significance of each multivariate term specified in the fitted 
manyglm model (Warton 2011). For each model term, it returns a test 
statistic as determined by the argument \code{test}, and a P-value calculated 
by resampling rows of the data using a method determined by the argument 
\code{resamp}. Of the four possible resampling methods, three (case, residual 
permutation and parametric boostrap) are described in more detail in Davison 
and Hinkley (1997, chapter 6), but the default (PIT-trap) is a new method (in 
review) which bootstraps probability integral transform residuals, and which 
we have found to give the most reliable Type I error rates. All methods 
involve resampling under the alternative hypothesis. These methods ensure 
approximately valid inference even when the mean-variance relationship or the 
correlation between variables has been misspecified. Standardized pearson 
residuals (see \code{\link{manyglm}} are currently used in residual 
permutation, and where necessary, resampled response values are truncated so 
that they fall in the required range (e.g. counts cannot be negative). 
However, this can introduce bias, especially for \code{family=binomial}, so 
we advise extreme caution using \code{perm.resid} for presence/absence data. 
If \code{resamp="none"}, p-values cannot be calculated, however the test 
statistics are returned. 

If you have a specific hypothesis of primary interest that you want to test, then you should use the \code{\link[mvabund:anova.manyglm]{anova.manyglm}} function, which can resample rows of the data under the null hypothesis and so usually achieves a better approximation to the true significance level.

For information on the different types of data that can be modelled using manyglm, see \code{\link{manyglm}}. To check model assumptions, use \code{\link{plot.manyglm}}.

Multivariate test statistics are constructed using one of three methods: a log-likelihood ratio statistic \code{test="LR"}, for example as in Warton et. al. (2012), or a Wald statistic \code{test="wald"} or a Score statistic \code{test="score"}. "LR" has good properties, but is only available when \code{cor.type="I"}.
%Also recommended, but only for multi-sample testing, is a Wald statistic using the variance-stabilising link function (e.g. \code{family=negative.binomial(varstab)}), which has been shown to have good power properties (Warton 2008b). However use of this link function is only reasonable for multi-sample tests, a special case where changing link function can be considered as a reparameterisation of the model (rather than a change of model).

The default Wald test statistic makes use of a generalised estimating equations (GEE) approach, estimating the covariance matrix of parameter estimates using a sandwich-type estimator that assumes the mean-variance relationship in the data is correctly specified and that there is an unknown but constant correlation across all observations. Such assumptions allow the test statistic to account for correlation between variables but to do so in a more efficient way than traditional GEE sandwich estimators (Warton 2008a). The common correlation matrix is estimated from standardized Pearson residuals, and the method specified by \code{cor.type} is used to adjust for high dimensionality.

The Wald statistic has problems for count data and presence-absence 
data when there are estimated means at zero (which usually means very large parameter estimates, check for this using \code{coef}). In such instances Wald statistics should not be used, Score or LR should do the job. 

%By default, \code{summary.manyglm} tests hypotheses concerning mean (absolute) abundance, but it is also possible to test hypotheses concerning community composition using the argument \code{compositional=TRUE}, when analysing count data, as in Warton and Guttorp (2007). This constructs Wald statistics to test hypotheses concerning relative abundance not absolute abundance, and the parameters that are returned in the fitted model are for the relative change in abundance of one variable as compared to another (as compared to the last abundance variable in the mvabund array).

The \code{summary.manyglm} function is designed specifically for high-dimensional data (that, is when the number of variables p is not small compared to the number of observations N). In such instances a correlation matrix is computationally intensive to estimate and is numerically unstable, so by default the test statistic is calculated assuming independence of variables (\code{cor.type="I"}). Note however that the resampling scheme used ensures that the P-values are approximately correct even when the independence assumption is not satisfied. However if it is computationally feasible for your dataset, it is recommended that you use \code{cor.type="shrink"} to account for correlation between variables, or \code{cor.type="R"} when p is small. The \code{cor.type="R"} option uses the unstructured correlation matrix (only possible when N>p), such that the standard classical multivariate test statistics are obtained. Note however that such statistics are typically numerically unstable and have low power when p is not small compared to N. 

The \code{cor.type="shrink"} option applies ridge regularisation (Warton (2008b)), shrinking the sample correlation matrix towards the identity, which improves its stability when p is not small compared to N. This provides a compromise between \code{"R"} and \code{"I"}, allowing us to account for correlation between variables, while using a numerically stable test statistic that has good properties. 

The shrinkage parameter is an attribute of the \code{\link{manyglm}} object. For a Wald test, the sample correlation matrix of the alternative model is used to calculate the test statistics. So \code{object$shrink.param} is used. For a Score test, the sample correlation matrix of the null model is used to calculate the test statistics. So \code{shrink.param} of the null model is used instead. If \code{cor.type=="shrink"} but \code{object$shrink.param} is not available, for example \code{object$cor.type!="shrink"}, then the shrinkage parameter will be estimated by cross-validation using the multivariate normal likelihood function (see \code{\link{ridgeParamEst}} and (Warton 2008b)) in the summary test.

Rather than stopping after testing for multivariate effects, it is often of interest to find out which response variables express significant effects. Univariate statistics are required to answer this question, and these are reported if requested. Setting \code{p.uni="unadjusted"} returns resampling-based univariate P-values for all effects as well as the multivariate P-values, whereas \code{p.uni="adjusted"} returns adjusted P-values (that have been adjusted for multiple testing), calculated using a step-down resampling algorithm as in Westfall & Young (1993, Algorithm 2.8). This method provides strong control of family-wise error rates, and makes use of resampling (using the method controlled by \code{resamp}) to ensure inferences take into account correlation between variables.
}
\value{
summary.manyglm returns an object of class "summary.manyglm", a list with components 
\item{call}{the component from \code{object}.}
\item{terms}{the terms object used.}
\item{family}{the component from \code{object}.}
\item{deviance}{the component from \code{object}.}
% \item{contrasts}{the component from \code{object}.}
\item{aic}{Akaike's \emph{An Information Criterion}, minus twice the 
maximized log-likelihood plus twice the number of coefficients 
(except for negative binomial and quasipoisson family, assuming 
that the dispersion is known).}
\item{df.residual}{the component from \code{object}.}
\item{null.deviance}{the component from \code{object}.}
\item{df.null}{the component from \code{object}.}
\item{devll}{minus twice the maximized log-likelihood}
\item{iter}{the number of iterations that were used in
\code{\link{manyglm}} for the estimation
of the model parameters.}
\item{p.uni}{the supplied argument.}
\item{nBoot}{the supplied argument.}
\item{resample}{the supplied argument.}
\item{na.action}{the na.action used in the \code{manyglm} object,
if applicable}
\item{show.residuals}{the supplied argument.}
\item{show.est}{the supplied argument.}
\item{compositional}{logical. Whether a test for 
  compositional effects was performed.}
% \item{show.ses}{the supplied argument.}
\item{test}{the supplied argument.}
\item{cor.type}{the supplied argument.}
\item{method}{the method used in \code{manyglm}. Either \code{"glm.fit"} or
\code{"manyglm.fit"}}
\item{theta.method}{the method used for the estimation of the nuisance parameter
theta.}
\item{manyglm.args}{a list of control parameters from \code{manyglm}.}
\item{rankX}{the rank of the design matrix.}
\item{covstat}{the supplied argument.}
\item{deviance.resid}{the deviance residuals.}
\item{est}{the estimated model coefficients}
\item{s.err}{the Scaled Variance}
\item{shrink.param}{
the shrinkage parameter. Either the value of the argument with the same name or
if this was not supplied the estimated shrinkage parameter.
}
\item{n.bootsdone}{the number of bootstrapping iterations that were done, i.e.
had no error.}
\item{coefficients}{
the matrix of coefficients, standard errors,
    z-values and p-values.  Aliased coefficients are omitted.
}
\item{stat.iter}{if the argument \code{stat.iter} is set to \code{TRUE} 
the test statistics in the resampling iterations.}
\item{statj.iter}{if the argument \code{stat.iter} is set to \code{TRUE}
the univariate test statistics in the resampling iterations.}
\item{aliased}{named logical vector showing if the original coefficients are aliased.} 
\item{dispersion}{either the supplied argument or the inferred/estimated
dispersion if the latter is \code{NULL}.
}
\item{df}{a 3-vector of the rank of the model and the number of residual degrees
of freedom, plus number of non-aliased coefficients.
}
\item{overall.n.bootsdone}{
the number of bootstrap iterations without errors that were done in the
overall test
}
\item{statistic}{a table containing test statistics, p values and degrees
of freedom for the overall test}
\item{overall.stat.iter}{if the argument \code{stat.iter} is set to \code{TRUE}
the test statistics of the overall tests in the resampling iterations.}
\item{overall.statj.iter}{if the argument \code{stat.iter} is set to \code{TRUE}
the univariate test statistics of the overall tests 
in the resampling iterations.}
\item{cov.unscaled}{
the unscaled (\code{dispersion = 1}) estimated covariance matrix of the
estimated coefficients.
}
\item{cov.scaled}{ditto, scaled by \code{dispersion}.}
\item{correlation}{
(only if the argument \code{show.cor = TRUE}.) The estimated correlations
of the estimated coefficients.
}
\item{symbolic.cor}{
(only if \code{show.cor = TRUE}.) The value of the argument \code{symbolic.cor}.
}
}
\references{
Warton D.I. (2011). Regularized sandwich estimators for analysis of high dimensional data using generalized estimating equations. \emph{Biometrics}, 67(1), 116-123.

Warton D.I. (2008a). Penalized normal likelihood and ridge regularization of 
correlation and covariance matrices. \emph{Journal of the American Statistical 
Association} 103, 340-349.

Warton D.I. (2008b). Which Wald statistic?  Choosing a parameterisation of the 
Wald statistic to maximise power in \emph{k}-sample generalised estimating 
equations. \emph{Journal of Statistical Planning and Inference}, 138, 3269-3282.

Warton D. I., Wright S., and Wang, Y. (2012). Distance-based multivariate analyses confound location and dispersion effects. \emph{Methods in Ecology and Evolution}, 3(1), 89-101.

Davison, A. C. and Hinkley, D. V. (1997) \emph{Bootstrap Methods and their Application},
Cambridge University Press, Cambridge.

Westfall, P. H. and Young, S. S. (1993) \emph{Resampling-based multiple testing.} 
John Wiley & Sons, New York.

Wu, C. F. J. (1986) Jackknife, Bootstrap and Other Resampling Methods in
Regression Analysis.
\emph{The Annals of Statistics} 14:4, 1261-1295.
}
% \note{...}
\author{
Yi Wang, David Warton <David.Warton@unsw.edu.au> and Ulrike Naumann.
}
\seealso{
\code{\link{manyglm}}, \code{\link[mvabund:anova.manyglm]{anova.manyglm}}.
}
\examples{
data(spider)
spiddat <- mvabund(spider$abund)

## Estimate the coefficients of a multivariate glm
glm.spid <- manyglm(spiddat[,1:3]~., data=spider$x, family="negative.binomial")

## Estimate the statistical significance of different multivariate terms in 
## the model, using the default settings of LR test, and 100 PIT-trap resamples
summary(glm.spid, show.time=TRUE) 

## Repeat with the parametric bootstrap and wald statistics 
summary(glm.spid, resamp="monte.carlo", test="wald", nBoot=300) 

}
\keyword{models}
\keyword{regression}
\keyword{multivariate}

