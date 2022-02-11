\name{manyany}
\alias{manyany}
\alias{print.manyany}
\title{Fitting Many Univariate Models to Multivariate Abundance Data}
\description{
  \code{manyany} is used to fit many univariate models (GLMs, GAMs, otherwise) to high-dimensional data, such as multivariate abundance data in ecology. This is the base model-fitting function - see \code{plot.manyany} for assumption checking, and \code{anova.manyany} for significance testing.
}
\usage{
manyany(formula, fn, family="negative.binomial", data, composition = FALSE, 
block = NULL, get.what="details", var.power=NA, na.action = "na.exclude", ...)
\method{print}{manyany}(x, digits = max(3L, getOption("digits") - 3L), ...) 
}
\arguments{
  \item{formula}{an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.  The details of model specification are given
    under Details.}
  \item{fn}{a character string giving the name of the function for the univariate
    model to be applied. e.g. "glm".}
  \item{family}{a description of the error distribution function to be used 
    in the model, either as a character string, a \code{\link{family}} object, or
    a list of such objects, one for each response variable in the dataset. Such a
    list enables the fitting of models with different distributions for different
    responses. See Details for the families currently supported.}
  \item{data}{an optional data frame containing predictor variables (a matrix is also acceptable).}
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
  \item{na.action}{Default set to \code{exclude} (for details see
    \code{\link{na.exclude}}) to avoid errors when NA's in predictors.}
  \item{\dots}{further arguments passed to the fitting function.}
  \item{x}{an object of class "manyany", usually, a result of a call to
  \code{\link[mvabund:manyany]{manyany}}.}
  \item{digits}{how many digits to include in the printed anova table.}
}
\details{

\code{manyany} can be used to fit the model type specified in \code{fn} to many variables
   simultaneously, a generalisation of \code{\link{manyglm}} to handle other model types.
   It should be able to handle any fixed effects modelling function that has 
   \code{predict} and \code{logLik} functions, and that accepts a \code{family} argument, 
   provided that the family is on our list (currently 'gaussian', 'poisson', 'binomial', 
   'negative.binomial' and 'tweedie', although models for ordinal data are also accepted 
   if using the \code{clm} function from the \code{ordinal} package). Models for \code{manyany}
   are specified symbolically, see for example the details section of \code{\link{lm}}
   and \code{\link{formula}}.

Unlike \code{manyglm}, this function accepts \code{family} functions as arguments
instead of just character strings, giving greater flexibility. For example, you could
use \code{family=binomial(link="probit")} to fit a model using the probit link, rather 
than being restricted to the default logit link or \code{cloglog} links available in \code{manyglm}.

A \code{data} argument is required, and it must be a dataframe.

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

spider = list(abund=spider$abund, x = as.matrix(spider$x))

## To fit a log-linear model assuming counts are negative binomial, via manyglm:
spidNB <- manyany(abund~x,"manyglm",data=spider,family="negative.binomial")

logLik(spidNB) # a number of generic functions are applible to manyany objects

## To fit a glm with complementary log-log link to presence/absence data:
PAdat = pmin(as.matrix(spider$abund),1) #constructing presence/absence dataset
spidPA <- manyany(PAdat~x,"glm",data=spider,family=binomial("cloglog"))
plot(spidPA)
# There are some wild values in there for the Pardmont variable (residuals >5 or <-8).
#The Pardmont model didn't converge, coefficients are a bit crazy:
coef(spidPA) 

# could try again using the glm2 package to fit the models, this fixes things up:
\dontrun{library(glm2)}
\dontrun{spidPA2<-manyany(PAdat~x,"glm",data=spider,family=binomial("cloglog"),method="glm.fit2")}
\dontrun{plot(spidPA2)} #looks much better.

## To simultaneously fit models to ordinal data using the ordinal package:
\dontrun{library(ordinal)}
## First construct an ordinal dataset:
\dontrun{spidOrd = spider$abund}
\dontrun{spidOrd[spider$abund>1 & spider$abund<=10]=2}
\dontrun{spidOrd[spider$abund>10]=3}
## Now fit a model using the clm function:
\dontrun{manyOrd=manyany(spidOrd~bare.sand+fallen.leaves,"clm",data=data.frame(spider$x))}
\dontrun{plot(manyOrd)}
}
\keyword{models}
\keyword{regression}
\keyword{multivariate}

