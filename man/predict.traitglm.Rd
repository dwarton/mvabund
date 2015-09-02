\name{predict.traitglm}
\alias{predict.traitglm}
\title{
Predictions from fourth corner model fits}
\description{
Obtains a prediction from a fitted fourth corner model object. 
}
 %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
\usage{
\method{predict}{traitglm}(object, newR=NULL, newQ=NULL, newL=NULL, type="response", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted object of class \code{traitglm}.
}
  \item{newR}{
A new data frame of environmental variables. If omitted, the original matrix of environmental variables is used.
}
  \item{newQ}{
A new data frame of traits for each response taxon. If omitted, the original matrix of traits is used.
}
  \item{newL}{
A new data frame of abundances (sites in rows, taxa in columns). This is only used if seeking predicted log-likelihoods. If omitted, the original abundances are used.
}
  \item{type}{
The type of prediction required. The default is predictions on the scale of the response variable, alternatives are \code{"logL"} for predictive log-likelihood, and "link" for linear predictors.
}
  \item{\dots}{
Further arguments passed to or from other methods.}
}
\details{
If \code{newR} and \code{newQ} are omitted, then as usual, predictions are based on the data used for the fit. Note that two types of predictions are possible in principle: predicting at new sites (by specifying a new set of environmental variables only, as \code{newR}) and predicting for new taxa (by specifying a new set of traits only, as \code{newQ}). Unfortunately, only predicting at new sites has been implemented at the moment! An issue with predicting to new taxa is that a main effect is included in the model for each taxon (by default), and the intercept would be unknown for a new species.

If predictive log-likelihoods are desired, a new data frame of abundances \code{newL} would need to be specified, whose rows correspond to those of \code{newR} and whose columns correspond to rows of \code{newQ}.
}
\value{
A matrix of predictions, with sites in rows and taxa in columns. 
}

\references{
Brown AM, Warton DI, Andrew NR, Binns M, Cassis G and Gibb H (2014) The fourth corner solution - using species traits to better understand how species traits interact with their environment, Methods in Ecology and Evolution 5, 344-352.
}
\author{
David I. Warton <David.Warton@unsw.edu.au>
}
\seealso{
\code{\link{traitglm}}}

 %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
\examples{
data(antTraits)

# fit a fourth corner model using negative binomial regression via manyglm:
ft=traitglm(antTraits$abund,antTraits$env,antTraits$traits,method="manyglm")
ft$fourth #print fourth corner terms

# predict to the first five sites
predict(ft, newR=antTraits$env[1:5,])

}
\keyword{models}
\keyword{regression}
