\name{predict.manyglm}
\alias{predict.manyglm}
\title{Predict Method for MANYGLM Fits}
\description{  
Obtains predictions and optionally estimates standard errors of those
predictions from a fitted manyglm object.
}
\usage{
% predict.manyglm
\method{predict}{manyglm}(object, newdata, type = c("link", "response",
    "terms"), se.fit = FALSE, dispersion = NULL, terms = NULL,
    na.action = na.pass, \dots)
}
\arguments{
  \item{object}{a fitted object of class inheriting from \code{"manyglm"}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict.  If omitted, the fitted linear predictors are used.}
  \item{type}{the type of prediction required.  The default is on the
    scale of the linear predictors; the alternative \code{"response"}
    is on the scale of the response variable.  Thus for a default
    binomial model the default predictions are of log-odds (probabilities
    on logit scale) and \code{type = "response"} gives the predicted
    probabilities.  The \code{"terms"} option returns a matrix giving the
    fitted values of each term in the model formula on the linear predictor
    scale.

    The value of this argument can be abbreviated.
  }
  \item{se.fit}{logical switch indicating if standard errors are required.}
  \item{dispersion}{the dispersion of the MANYGLM fit to be assumed in
    computing the standard errors.  If omitted, that returned by
    \code{summary} applied to the object is used.}
  \item{terms}{with \code{type="terms"} by default all terms are returned.
    A character vector specifies which terms are to be returned}
  \item{na.action}{function determining what should be done with missing
    values in \code{newdata}. The default is to predict \code{NA}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
predict.manyglm refits the model using glm before making predictions.  In rare (usually pathological) cases this may lead to differences in predictions as compared to what would be expected if using the manyglm coefficients directly.

If \code{newdata} is omitted the predictions are based on the data
  used for the fit.  In that case how cases with missing values in the
  original fit is determined by the \code{na.action} argument of that
  fit.  If \code{na.action = na.omit} omitted cases will not appear in
  the residuals, whereas if \code{na.action = na.exclude} they will
  appear (in predictions and standard errors), with residual value
  \code{NA}.  See also \code{\link{napredict}}.
}
\value{
  If \code{se = FALSE}, a matrix of predictions or an array of 
  predictions and bounds.  
  If \code{se = TRUE}, a list with components
  \item{fit}{the predictions}
  \item{se.fit}{estimated standard errors}
  \item{residual.scale}{a scalar giving the square root of the
    dispersion used in computing the standard errors.}
}
% \references{...}
% \note{...}
\author{Ulrike Naumann, Yi Wang and David Warton <David.Warton@unsw.edu.au>.}
\seealso{
\code{\link{manyglm}}.
}
\examples{
data(spider)
spiddat <- mvabund(spider$abund)
Y <- spiddat[1:20,]
X <- spider$x[1:20,]
glm.spid.poiss <- manyglm(Y~soil.dry+bare.sand, family="poisson", data=X)
glm.spid.poiss$data = X
newdata <- spider$x[21:28,]
predict(glm.spid.poiss, newdata)
pred.w.plim <- predict(glm.spid.poiss, newdata, interval="prediction")
pred.w.clim <- predict(glm.spid.poiss, newdata, interval="confidence")
}
\keyword{models}
\keyword{multivariate}
\keyword{regression}

