\name{manylm}
\alias{manylm}
\title{Fitting Linear Models for Multivariate Abundance Data}
\description{
\code{manylm} is used to fit multivariate linear models
to high-dimensional data, such as multivariate abundance data in ecology.

This is the base model-fitting function - see \code{plot.manylm} for 
assumption checking, and \code{anova.manylm} or \code{summary.manylm} 
for significance testing.
}
\usage{
manylm(
   formula, data=NULL,  subset=NULL, weights=NULL, 
   na.action=options("na.action"),  method="qr", model=FALSE, 
   x=TRUE, y=TRUE, qr=TRUE, singular.ok=TRUE, contrasts=NULL, 
   offset, test="LR" , cor.type= "I", shrink.param=NULL, 
   tol=1.0e-5, \dots)
}
\arguments{
  \item{formula}{
    an object of class \code{"formula"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted. The details of model specification are given
    under Details.
  }

  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{as.data.frame} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{manylm} is called.}

  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}

  \item{weights}{an optional vector of weights to be used in the fitting
    process.  Should be \code{NULL} or a numeric vector.
    If non-null, weighted least squares is used with weights
    \code{weights} (that is, minimizing \code{sum(weights*e^2)}); otherwise
    ordinary least squares is used.}

  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{options}, and is
    \code{na.fail} if that is unset.  The \sQuote{factory-fresh}
    default is \code{na.omit}.  Another possible value is
    \code{NULL}, no action.  Value \code{na.exclude} can be useful.}

  \item{method}{the method to be used; for fitting, currently only
    \code{method = "qr"} is supported; \code{method = "model.frame"} returns
    the model frame (the same as with \code{model = TRUE}, see below).}

  \item{model, x, y, qr}{logicals.  If \code{TRUE} the corresponding
    components of the fit (the model frame, the model matrix, the
    response, the QR decomposition) are returned. }

  \item{singular.ok}{logical. If \code{FALSE} (the default in S but
    not in \R) a singular fit is an error.}

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}

  \item{offset}{this can be used to specify an \emph{a priori}
    known component to be included in the linear predictor
    during fitting.  This should be \code{NULL} or a numeric vector of
    length either one or equal to the number of cases.
    One or more \code{offset} terms can be included in the
    formula instead or as well, and if both are specified their sum is
    used.  See \code{model.offset}.}

  \item{test}{choice of test statistic. Can be one of "LR" (default) = likelihood ratio statistic "F" = Lawley-Hotelling trace statistic \code{NULL} = no test
This parameter is merely stored in \code{manylm}, and will be used as the default value of \code{test} in subsequent functions for inference.}

\item{cor.type}{structure imposed on the estimated correlation matrix under the fitted model. Can be "I"(default), "shrink", or "R". See anova.manylm for details of its usage. This parameter will be used as the default value of \code{cor.type} in subsequent functions for inference.}

\item{shrink.param}{shrinkage parameter to be used if \code{cor.type="shrink"}. This parameter will be used as the default value of \code{shrink.param} in subsequent functions for inference.}

\item{tol}{the tolerance used in estimations.}
\item{\dots}{additional arguments to be passed to the low level regression fitting functions (see below).}
}

\details{
  Models for \code{manylm} are specified symbolically. For details on this 
  compare the details section of \code{lm} and \code{formula}. If the formula 
  includes an \code{offset}, this is evaluated and subtracted from the 
  response.   \cr
  See \code{model.matrix} for some further details. The terms in
  the formula will be re-ordered so that main effects come first,
  followed by the interactions, all second-order, all third-order and so
  on: to avoid this pass a \code{terms} object as the formula (see
  \code{aov} and \code{demo(glm.vr)} for an example). \cr
  A formula has an implied intercept term.  To remove this use either
  \code{y ~ x - 1} or \code{y ~ 0 + x}.  See \code{formula} for
  more details of allowed formulae. \cr
  \code{manylm} calls the lower level function \code{\link{manylm.fit}} 
  or \code{\link{manylm.wfit}} for the actual numerical computations. 
  For programming only, you may consider doing likewise. \cr
  All of \code{weights}, \code{subset} and \code{offset} are evaluated
  in the same way as variables in \code{formula}, that is first in
  \code{data} and then in the environment of \code{formula}.\cr
  For details on arguments related to hypothesis testing (such as \code{cor.type}
  and \code{resample}) see \code{\link{summary.manylm}} or
  \code{\link{anova.manylm}}.
}
\value{
 \code{manylm} returns an object of \code{c("manylm", "mlm", "lm")} for multivariate
 formula response and of of class \code{c("lm")} for univariate response.

A \code{manylm} object is a list containing at least the following components:
  % \item{assign}{the assign attribute of the design matrix}
  \item{coefficients}{a named matrix of coefficients}
  \item{residuals}{the residuals matrix, that is response minus fitted values.}
  \item{fitted.values}{the matrix of the fitted mean values.}
  \item{rank}{the numeric rank of the fitted linear model.}
  %  \item{qr}{}
  \item{weights}{(only for weighted fits) the specified weights.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{hat.X}{the hat matrix.}
  \item{txX}{the matrix \code{(t(x)\%*\%x)}.}
  \item{test}{the \code{test} argument supplied.}
  \item{cor.type}{the \code{cor.type} argument supplied.}
  \item{resample}{the \code{resample} argument supplied.}
  \item{nBoot}{the \code{nBoot} argument supplied.}
  \item{call}{the matched call.}
  \item{terms}{the \code{terms} object used.}
  % \item{contrasts}{(only where relevant) the contrasts used.}
  \item{xlevels}{(only where relevant) a record of the levels of the
    factors used in fitting.}
  \item{model}{if requested (the default), the model frame used.}
  \item{offset}{the offset used (missing if none were used).}
  \item{y}{if requested, the response matrix used.}
  \item{x}{if requested, the model matrix used.}

  In addition, non-null fits will have components \code{assign} and 
  (unless not requested) \code{qr} relating to the linear
  fit, for use by extractor functions such as \code{summary.manylm}.
}
% \references{...}
% \note{...}
\author{
Yi Wang, Ulrike Naumann and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{anova.manylm}}, \code{\link{summary.manylm}}, \code{\link{plot.manylm}}}
\examples{
data(spider)
spiddat <- log(spider$abund+1)
lm.spider <- manylm(spiddat~.,data=spider$x)
lm.spider

#Then use the plot function for diagnostic plots, and use anova or summary to
#evaluate significance of different model terms.
}
\keyword{models}
\keyword{multivariate}
\keyword{regression}

