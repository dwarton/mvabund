\name{anova.traitglm}
\alias{anova.traitglm}
\title{Testing for a environment-by-trait (fourth corner) interaction by analysis of deviance}
\description{Returns an analysis of deviance from a fourth corner model, as computed using \code{\link{traitglm}}, typically to test for an environment-by-trait interaction. Slowly! This function works via \code{\link{anova.manyglm}}, which uses row-resampling for inference, and it only applies to \code{\link{traitglm}} objects that have been fitted using the (default) \code{\link{manyglm}} function.}

\usage{
 %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
\method{anova}{traitglm}(object, \dots, nBoot=99, resamp="pit.trap", test="LR",
                             block = NULL, show.time="all", bootID=NULL)
}
\arguments{
  \item{object}{
A fitted object of class \code{traitglm} and class \code{manyglm}. 
}
  \item{\dots}{
Additional \code{traitglm} objects, fitted using the \code{formula} argument.
}
  \item{nBoot}{
The number of bootstrap iterations. Default is 99 (NOTE: you should increase this for later runs!)
}
  \item{resamp,test,bootID}{
Arguments as in \code{\link{anova.manyglm}}, to control resampling method (\code{resamp}), test statistic (\code{test}) and whether or not a matrix of bootstrap resamples is manually entered (\code{bootID}).
}
\item{block}{
A factor specifying the sampling level to be resampled. Default is resampling sites (which still involves passing a blocking variable to \code{anova.manyglm}, to keep all rows of the original abundance data together in resamples).
}
\item{show.time}{Whether to display timing information for the resampling procedure: this is advisable, as resampling fourth corner models many times can take a while. The default value \code{"all"} shows all timing information, \code{"total"} shows only the overall time taken for the tests, and \code{"none"} shows none. }
}

\details{
There are two possible uses of this function, depending whether one \code{traitglm} object is specified or multiple objects.

If one \code{traitglm} object is specified, the \code{anova.traitglm} function returns a table summarising the statistical significance of the fourth corner terms in a model, that is, the interaction between environment and traits in predicting abundance across taxa and sites. All environment-by-trait interaction terms from the model are simultaneously tested.

If two or more nested \code{traitglm} objects are specified, and each has been fitted using a \code{formula} argument to the same set of datasets, then sequential test statistics (and P values) are returned for each additional model fit.

All \code{\link{traitglm}} models must be fitted using the \code{\link{manyglm}} function, which is its default behaviour, in order to access the \code{\link{anova.manyglm}}, which does most of the work. See \code{\link{anova.manyglm}} for details on how resampling is done, and options for arguments controlling the test statistic (via \code{test}) and the resampling method (via \code{resamp}). Because \code{\link{traitglm}} models are fitted by first vectorising the data into a univariate model, arguments such as \code{p.uni} and \code{cor.type} are redundant.

\code{\link{traitglm}} fits a single model to abundances across all sites and taxa at the same time, meaning the vector of abundances is typically pretty long, and the design matrix explaining how abundance varies across taxa and sites is typically pretty large. So resampling can take yonks. Hence the default number of resamples has been set at \code{nBoot=99}, but please consider increasing this once you have a sense for how long it will take to run (scales roughly linearly with \code{nBoot}).
}
\value{
A list of values as returned by \code{\link{anova.manyglm}}, of which the most relevant element is \code{table} (the analysis of deviance table).
}

\section{Warning}{
The comparison between two or more models by \code{anova.traitglm} will only be valid if they are fitted to the same dataset. This may be a problem if there are missing values and R's default of \code{na.action = na.omit} is used. 
}

\references{
Warton DI, Shipley B & Hastie T (2015) CATS regression - a model-based approach to studying trait-based community assembly, Methods in Ecology and Evolution 6, 389-398.
}
\author{
David I. Warton <David.Warton@unsw.edu.au>
}
\seealso{
\code{\link[mvabund:anova.manyglm]{anova.manyglm}},\code{\link{traitglm}}
}
 %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
\examples{
data(antTraits)

# we'll fit a small fourth corner model, to a subset of the antTraits data.
# first select only species present in at least 25% of plots:
abSum = apply(antTraits$abund>0,2,mean)
ab = antTraits$abund[,abSum>0.25]
tr = antTraits$traits[abSum>0.25,]

# now fit the fourth corner model, only as a function of a couple of traits and env variables:
ft=traitglm(ab,antTraits$env[,1:3],data.frame(tr$Weber,tr$Femur))
anova(ft,nBoot=39)
# Note you should refit with more bootstrap samples (e.g. 999), should take <2 minutes to run

# for an example using anova.traitglm for multiple fits, uncomment the following lines:
# ft2=traitglm(antTraits$abund,antTraits$env[,3:4],antTraits$traits[,c(1,3)],
#   formula=~1,composition=TRUE) #no fourth corner terms
# ft3=traitglm(antTraits$abund,antTraits$env[,3:4],antTraits$traits[,c(1,3)],
#   formula=~Shrub.cover:Femur.length+Shrub.cover:Pilosity,composition=TRUE) #shrub interactions
# ft4=traitglm(antTraits$abund,antTraits$env[,3:4],antTraits$traits[,c(1,3)],
#   formula=~Shrub.cover:Femur.length+Shrub.cover:Pilosity+Volume.lying.CWD:Femur.length+
#   Volume.lying.CWD:Pilosity, composition=TRUE) #all interactions only
# anova(ft2,ft3,ft4) # Shrub interactions not significant but CWD interactions are.
}
\keyword{models}
\keyword{regression}
\keyword{multivariate}

