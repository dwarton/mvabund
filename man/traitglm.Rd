\name{traitglm}
\alias{traitglm}
\title{
Fits a fourth corner model for abundance as a function of environmental variables and species traits.
}
\description{
Fits a fourth corner model - a model to study how variation in environmental response across taxa can be explained by their traits. The function to use for fitting can be (pretty well) any predictive model, default is a generalised linear model, another good option is to add a LASSO penalty via \code{glm1path}. Can handle overdispersed counts via  \code{family="negative.binomial"}, which is the default \code{family} argument.
}
\usage{
 %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
traitglm(L, R, Q = NULL, family="negative.binomial", formula = NULL, method = "manyglm",
            composition = FALSE, col.intercepts = TRUE, ...)
}
\arguments{
  \item{L}{
A data frame (or matrix) containing the abundances for each taxon (columns) across all sites (rows).
}
  \item{R}{
A data frame (or matrix) of environmental variables (columns) across all sites (rows).
}
  \item{Q}{
A data frame (or matrix) of traits (columns) across all taxa (rows). If not specified, a different environmental response will be specified for each taxon.
}
  \item{family}{
The family of the response variable, see \code{\link{family}}. Negative binomial with unknown overdispersion and a log-link can be specified as "negative.binomial", and is the default.
}
  \item{formula}{
A one-sided formula specifying exactly how to model abundance as a function of environmental and trait variables (as found in \code{R} and \code{Q} respectively).  Default is to include all terms additively, with quadratics for quantitative terms, and all environment-by-trait interactions.
}
  \item{method}{
The function to use to fit the model. Default is \code{\link{manyglm}}, some other available options are \code{\link{glm1path}}, \code{\link{cv.glm1path}} for LASSO-penalised fits, but in principle any model-fitting function that accepts formula input and a family argument should work.
}
\item{composition}{logical. TRUE includes a row effect in the model, adjusting for different sampling intensities across different samples. This can be understood as a compositional term in the sense that all other terms then model relative abundance at a site. FALSE (default) does not include a row effect, hence the model is of absolute abundance.
}
\item{col.intercepts}{logical. TRUE (default) includes a column effect in the model, to adjust for different levels of abundance of different response (column) variables. FALSE removes this column effect.
}
%  \item{fourthPlot}{
%Logical, default (\code{TRUE}) will produce a lattice plot of fourth corner coefficients, that is, of interactions between environmental and %trait variables in predicting abundance. Otherwise nothing is plotted.
%}
  \item{\dots}{
Arguments passed to the function specified at \code{method} that will be used to fit the model.}
}
\details{
This function fits a fourth corner model, that is, a model to predict abundance across several taxa (stored in \code{L}) as a function of environmental variables (\code{R}) and traits (\code{Q}). The environment-trait interaction can be understood as the fourth corner, giving the set of coefficients that describe how environmental response across taxa varies as traits vary. A species effect is include in the model (i.e. a different intercept term for each species), so that traits are used to explain patterns in relative abundance across taxa not patterns in absolute abundance.

The actual function used to fit the model is determined by the user through the \code{method} argument. The default is to use \code{\link{manyglm}} to fit a GLM, although for predictive modelling, it might be better to use a LASSO penalty as in \code{\link{glm1path}} and \code{\link{cv.glm1path}}. In \code{\link{glm1path}}, the penalty used for BIC calculation is \code{log(dim(L)[1])}, i.e. log(number of sites).

The model is fitted by vectorising \code{L} then constructing a big matrix from repeated values of \code{R}, \code{Q}, their quadratic terms (if required) and interactions. Hence this function will hit memory issues if any of these matrices are large, and can slow down (especially if using \code{\link{cv.glm1path}}). If \code{formula} is left unspecified, the design matrix is constructed using all environmental variables and traits specified in \code{R} and \code{Q}, and quadratic terms for any of these variables that are quantitative, and all environment-trait interactions, after standardising these variables. Specifying a one-sided \code{formula} as a function of the variables in \code{R} and \code{Q} would instead give the user control over the precise model that is fitted, and drops the internal standardisations. The arguments \code{composition} and \code{col.intercepts} optionally add terms to the model for row and column total abundance, irrespective of whether a \code{formula} has been specified.

Note: when specifying a \code{formula}, if there are no penalties on coefficients (as for \code{\link{manyglm}}), then main effects for \code{R} can be excluded if including row effects (via \code{composition=TRUE}), and main effects for \code{Q} can be excluded if including column effects (via \code{col.intercepts=TRUE}), because those terms are redundant (trying to explain main effects for row/column when these main effects are already in the model). If using penalised likelihood (as in \code{\link{glm1path}} and \code{\link{cv.glm1path}}) or a random effects model, by all means include main effects as well as row/column effects, and the penalties will sort out which terms to use.

If trait matrix \code{Q} is not specified, default behaviour will fit a different environmental response for each taxon (and the outcome will be very similar to \code{manyglm(L~R)}). This can be understood as a fourth corner model where species identities are used as the species traits (i.e. no attempt is made to explain differences across species).

These functions inherit default behaviour from their fitting functions. e.g. use \code{plot} for a Dunn-Smyth residual plot from a traits model fitted using \code{\link{manyglm}} or \code{\link{glm1path}}.
}
\value{Returns a \code{traitglm} object, a list that contains at least the following components:
\describe{
\item{\dots}{Exactly what is included in output depends on the fitting function - by default, a \code{\link{manyglm}} object is returned, so all usual \code{manyglm} output is included (coefficients, residuals, deviance, etc).}
\item{family}{A \code{family} object matching the final model.}
\item{fourth.corner}{A matrix of fourth corner coefficients. If \code{formula} has been manually entered, this will be a vector not a matrix.}
\item{R.des}{The reduced-size design matrix for environmental variables, including further arguments:
\describe{
\item{X}{Data frame of (possibly standardised) environmental variables}
\item{X.squ}{A data frame containing the leading term in a quadratic expression (where appropriate) for environmental variables}
\item{var.type}{A vector with the same dimension as the number of columns of X, listing the type of ecah enviromental variable (\code{"quantitative"}" or \code{"factor"}")}
\item{coefs}{Coefficients used in transforming variables to orthogonality. These are used later to make predictions.}
}
}
\item{Q.des}{The reduced-size design matrix for traits, set up as for \code{R.des}.}
\item{spp.penalty}{For LASSO fits: a vector of the same length as the final design matrix, indicating which variables had a penalty imposed on them in model fitting.}
\item{L}{The data frame of abundances specified as input.}
\item{any.penalty}{Logical, is any penalty applied to parameters at all (not if using a \code{manyglm} fit).}
\item{scaling}{A list of coefficients describing the standaridsations of variables used in analyses. Stored for use later if making predictions.}
\item{call}{The original call \code{traitglm} call.}
}
}

\references{
Brown AM, Warton DI, Andrew NR, Binns M, Cassis G and Gibb H (2014) The fourth corner solution - using species traits to better understand how species traits interact with their environment, Methods in Ecology and Evolution 5, 344-352.

Warton DI, Shipley B & Hastie T (2015) CATS regression - a model-based approach to studying trait-based community assembly, Methods in Ecology and Evolution 6, 389-398.
}
\author{
David I. Warton <David.Warton@unsw.edu.au>
}
\seealso{
\code{\link{glm1path}}, \code{\link{glm1}}, \code{\link{manyglm}}, \code{\link{family}}, \code{\link[mvabund:residuals.manyglm]{residuals.manyglm}}, \code{\link[mvabund:plot.manyany]{plot.manyany}}
}

 %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
\examples{
data(antTraits)

ft=traitglm(antTraits$abund,antTraits$env,antTraits$traits,method="manyglm")
ft$fourth #print fourth corner terms

# for a pretty picture of fourth corner coefficients, uncomment the following lines:
# library(lattice)
# a        = max( abs(ft$fourth.corner) )
# colort   = colorRampPalette(c("blue","white","red")) 
# plot.4th = levelplot(t(as.matrix(ft$fourth.corner)), xlab="Environmental Variables",
#   ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
#   scales = list( x= list(rot = 45)))
# print(plot.4th)

plot(ft) # for a Dunn-smyth residual plot
qqnorm(residuals(ft)); abline(c(0,1),col="red") # for a normal quantile plot.

# predict to the first five sites
predict(ft,newR=antTraits$env[1:5,])

# refit using LASSO and less variables, including row effects and only two interaction terms:
ft1=traitglm(antTraits$abund,antTraits$env[,3:4],antTraits$traits[,c(1,3)],
      formula=~Shrub.cover:Femur.length+Shrub.cover:Pilosity,composition=TRUE,method="glm1path")
ft1$fourth #notice LASSO penalty has one interaction to zero

}
\keyword{models}
\keyword{regression}
