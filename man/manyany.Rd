\name{manyany}
\alias{manyany}
\title{Fitting Many Univariate Models to Multivariate Abundance Data}
\description{
  \code{manyany} is used to fit many univariate models (GLMs, GAMs, otherwise) to high-dimensional data, such as multivariate abundance data in ecology. This is the base model-fitting function - see \code{plot.manyany} for assumption checking, and \code{anova.manyany} for significance testing.
}
\usage{
manyany(fn, yMat, formula, data, family="negative.binomial", composition = FALSE, 
block = NULL, get.what="details", var.power=NA, ...)
}
\arguments{
  \item{fn}{a character string giving the name of the function for the univariate
    model to be applied. e.g. "glm".}
  \item{yMat}{a matrix of response variables, e.g. multivariate abundances.}
  \item{formula}{an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.  The details of model specification are given
    under Details.}
  \item{data}{a data frame containing predictor variables (a matrix is also acceptable).
    This is REQUIRED and needs to have more than one variable in it (even if only one
    is used in the model).}
  \item{family}{a description of the error distribution function to be used 
    in the model, either as a character string, a \code{\link{family}} object, or
    a list of such objects, one for each response variable in the dataset. Such a
    list enables the fitting of models with different distributions for different
    responses. See Details for the families currently supported.}
  \item{composition}{logical. FALSE (default) fits a separate model to each species.
    TRUE fits a single model to all variables, including site as a row effect, such
    that all other terms model relative abundance (compositional effects).  }
  \item{block}{a factor specifying the sampling level to be resampled. Default is
    resampling rows (if composition=TRUE in the manyany command, this means resampling
    rows of data as originally sent to manyany).}
  \item{get.what}{what to return from each model fit: "details" (default) includes
    predicted values and residuals in output, "models" also returns the fitted objects
    for each model, "none" returns just the log-likelihood (mostly for internal use).}
  \item{var.power}{the power parameter, if using the tweedie distribution.}
  \item{\dots}{further arguments passed to the fitting function.}
}
\details{

\code{manyany} can be used to fit the specified model type to many variables
   simultaneously, a generalisation of \code{\link{manyglm}}. It should be able to handle
   any fixed effects modelling function that has \code{predict} and \code{logLik}
   functions, and that accepts a \code{family} argument, provided that the family
   is on our list (currently 'gaussian', 'poisson', 'binomial', 'negative.binomial'
   and 'tweedie', although models for ordinal data are also accepted if using the
   \code{clm} function form the \code{ordinal} package). Models for \code{manyany}
   are specified symbolically, see for example the details section of \code{\link{lm}}
   and \code{\link{formula}}.

Unlike \code{manyglm}, this function accepts \code{family} functions as arguments
instead of just character strings, giving greater flexibility. For example, you can
use \code{family=binomial(link="cloglog")} to fit a model using the complementary
log-log link, rather than being restricted to the default logit link.

A \code{data} argument is required, and it must be a dataframe containing more than
one object.  It need not contain that matrix of response variables, that is specified
separately as \code{yMat}.

Setting \code{composition=TRUE} enables compositional analyses, where predictors are
used to model relative abundance rather than mean abundance. This is achieved by
vectorising the response matrix and fitting a single model across all variables, with
a row effect to account for differences in relative abundance across rows.
The default \code{composition=FALSE} just fits a separate model for each variable.
}

\value{
  \code{manyany} returns an object inheriting from \code{"manyany"}.

  The function \code{anova} (i.e. \code{\link[mvabund:anova.manyany]{anova.manyany}}) will produce a significance test comparing two \code{manyany} objects.
  Currently there is no \code{summary} resampling function for objects of this class.

  The generic accessor functions \code{fitted.values}, \code{residuals}, \code{logLik}, \code{AIC}, \code{plot}
  can be used to extract various useful features of the value returned by \code{manyany}.

  An object of class \code{"manyany"} is a list containing at least the
  following components:

  \item{logL}{a vector of log-likelihood terms for each response variable in the fitted model.}
  \item{fitted.values}{the matrix of fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
  \item{residuals}{the matrix of probability integral transform (PIT) residuals.
     If the fitted model is a good fit, these will be approximately standard uniformly distributed.}
  \item{linear.predictor}{the linear fit on link scale. But for ordinal models fitted using \code{clm}, these values are for the first category only.}
  \item{family}{a vector of \code{\link{family}} arguments, one for each response variable.}
  \item{call}{the matched call.}
  \item{model}{the \code{model.frame} from the model for the last response variable.}
  \item{terms}{a list of \code{terms} from the model for the last response variable.}
}

\references{

Warton D. I., Wright S., and Wang, Y. (2012). Distance-based multivariate analyses confound location and dispersion effects. \emph{Methods in Ecology and Evolution}, 3(1), 89-101.

}

\author{
David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link[mvabund:anova.manyglm]{anova.manyany}}, \code{\link[mvabund:residuals.manyany]{residuals.manyany}}, \code{\link[mvabund:plot.manyany]{plot.manyany}}.
}
\examples{
data(spider)
abund <- spider$abund
X <- as.matrix(spider$x)

## To fit a log-linear model assuming counts are negative binomial, via manyglm:
spidNB <- manyany("manyglm",abund,data=X,abund~X,family="negative.binomial")

logLik(spidNB) # a number of generic functions are applible to manyany objects

## To fit a glm with complementary log-log link to presence/absence data:
PAdat = pmin(as.matrix(abund),1) #constructing presence/absence dataset
spidPA <- manyany("glm",PAdat,data=X,PAdat~X,family=binomial("cloglog"))
plot(spidPA)
# There are some wild values in there for the Pardmont variable (residuals >5 or <-8).
#The Pardmont model didn't converge, coefficients are a bit crazy:
coef(spidPA) 

# Can try again using the glm2 package to fit the models, this fixes things up:
# library(glm2)
# spidPA2 <- manyany("glm",PAdat,data=X,PAdat~X,family=binomial("cloglog"),method="glm.fit2")
# plot(spidPA2) #looks much better.

## To simultaneously fit models to ordinal data using the ordinal package:
# library(ordinal)
## First construct an ordinal dataset:
# spidOrd = abund
# spidOrd[abund>1 & abund<=10]=2
# spidOrd[abund>10]=3
# for(iVar in 1:dim(spidOrd)[2])
#   spidOrd[,iVar]=factor(spidOrd[,iVar])
##Now fit a model using the clm function:
# manyOrd=manyany("clm",spidOrd,abund~bare.sand+fallen.leaves,data=X)
# plot(manyOrd)

}
\keyword{models}
\keyword{regression}
\keyword{multivariate}

