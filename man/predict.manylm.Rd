\name{predict.manylm}
% add details and note
\alias{predict.manylm}
\title{Model Predictions for Multivariate Linear Models}
\description{\code{predict.manylm} is a function for predictions from the
  result of the model fitting function \code{manylm}.  }
\usage{
% predict.manylm
\method{predict}{manylm}(object, newdata=NULL, se.fit = FALSE, 
    type = c("response", "terms"), terms = NULL, na.action = na.pass, ...)
}
\arguments{
\item{object}{object of class inheriting from \code{manylm}.}
\item{newdata}{an optional data frame in which to look for variables with which
to predict. If omitted, the fitted values are used.}
\item{se.fit}{a switch indicating if standard errors are required.}
\item{type}{type of prediction (response or model term), Possible values:
"response", "terms".}
\item{terms}{if type="terms", which terms (default is all terms).}
\item{na.action}{function determining what should be done with missing values in
newdata. The default is to predict \code{NA}.}
\item{\dots}{further arguments passed to or from other methods.}
}
\details{ \code{predict.manylm} produces predicted values, obtained by evaluating
  the regression function in the frame \code{newdata} (which defaults to
  \code{model.frame(object)}.  If the logical \code{se.fit} is
  \code{TRUE}, standard errors of the predictions are calculated.  If
  the numeric argument \code{scale} is set (with optional \code{df}), it
  is used as the residual standard deviation in the computation of the
  standard errors, otherwise this is extracted from the model fit.
  Setting \code{intervals} specifies computation of confidence or
  prediction (tolerance) intervals at the specified \code{level}, sometimes 
  referred to as narrow vs. wide intervals.

  If the fit is rank-deficient, some of the columns of the design matrix
  will have been dropped.  Prediction from such a fit only makes sense
  if \code{newdata} is contained in the same subspace as the original
  data.  That cannot be checked accurately, so a warning is issued.

  If \code{newdata} is omitted the predictions are based on the data
  used for the fit.  In that case how cases with missing values in the
  original fit is determined by the \code{na.action} argument of that
  fit.  If \code{na.action = na.omit} omitted cases will not appear in
  the residuals, whereas if \code{na.action = na.exclude} they will
  appear (in predictions, standard errors or interval limits),
  with residual value \code{NA}.  See also \code{\link{napredict}}.

  The prediction intervals are for a single observation at each case in
  \code{newdata} (or by default, the data used for the fit) with error
  variance(s) \code{pred.var}. This can be a multiple of \code{res.var},
  the estimated
  value of \eqn{\sigma^2}{sigma^2}: the default is to assume that future
  observations have the same error variance as those
  used for fitting. If \code{weights} is supplied, the inverse of this
  is used as a scale factor. For a weighted fit, if the prediction
  is for the original data frame, \code{weights} defaults to the weights
  used for the  model fit, with a warning since it might not be the
  intended result. If the fit was weighted and newdata is given, the
  default is to assume constant prediction variance, with a warning.
}
\value{
\code{predict.manylm} produces a matrix of predictions or if \code{interval}
is set an array of predictions and bounds, where the first dimension has the names:
\code{fit}, \code{lwr}, and \code{upr}.
If \code{se.fit} is \code{TRUE}, a list with the following components is returned:
\item{fit}{vector or matrix as above}
\item{se.fit}{a matrix with the standard errors of the predicted means}
\item{residual.scale}{vector or matrix as a vector of residual standard deviations}
\item{df}{numeric, the degrees of freedom for residual}
}
\note{
  Variables are first looked for in \code{newdata} and then searched for
  in the usual way (which will include the environment of the formula
  used in the fit).  A warning will be given if the
  variables found are not of the same length as those in \code{newdata}
  if it was supplied.
  
  Offsets specified by \code{offset} in the fit by \code{\link{lm}}
  will not be included in predictions, whereas those specified by an
  offset term in the formula will be.

  Notice that prediction variances and prediction intervals always refer
  to \emph{future} observations, possibly corresponding to the same
  predictors as used for the fit. The variance of the \emph{residuals}
  will be smaller.
  
  Strictly speaking, the formula used for prediction limits assumes that
  the degrees of freedom for the fit are the same as those for the
  residual variance.  This may not be the case if \code{res.var} is
  not obtained from the fit. 
}
% \references{...}
% \note{...}
\author{Ulrike Naumann, Yi Wang and David Warton <David.Warton@unsw.edu.au>.}
\seealso{
\code{\link{manylm}}.
}
\examples{
data(spider)
spiddat <- mvabund(spider$abund[1:20, ])
dat = data.frame(spider$x[1:20, ])
manylm.fit <- manylm(spiddat~soil.dry+bare.sand, data=dat)
predict(manylm.fit)
predict(manylm.fit, se.fit = TRUE)

new <- data.frame(spider$x[21:28, ])
predict(manylm.fit, new, se.fit = TRUE)

}
\keyword{methods}
\keyword{multivariate}
