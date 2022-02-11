\name{mvabund-package}

\alias{mvabund-package}

\title{Statistical methods for analysing multivariate abundance data }

\description{
This package provides tools for a model-based approach to the analysis of multivariate abundance data in ecology (Warton 2011), where 'abundance' should be interpreted loosely - as well as counts you could have presence/absence, ordinal or biomass (via \code{\link{manyany}}), etc.

There are graphical methods for exploring the properties of data and the community-environment association, flexible regression methods for estimating and making robust inferences about the community-environment association, 'fourth corner models' to explain environmental response as a function of traits, and diagnostic plots to check the appropriateness of a fitted model (Wang et. al 2012).

There is an emphasis on design-based inferences about these models, e.g. bootstrapping rows of residuals via \code{anova} calls, or cross-validation across rows, to make multivariate inferences that are robust to failure of assumptions about correlation. Another emphasis is on presenting diagnostic tools to check assumptions, especially via residual plotting.
}

\details{
The key functions available in this package are the following.

\bold{For graphical display of the data:}
    \describe{
      \item{\code{\link{plot.mvabund}}}{ draw a range of plots for Multivariate Abundance Data }
      \item{\code{\link{boxplot.mvabund}}}{ draw a range of plots of Model Formulae for Multivariate Abundance Data }
      \item{\code{\link{meanvar.plot}}}{ draw mean-variance plots for Multivariate Abundance Data }
    }

\bold{For estimating and displaying Linear Models:}
    \describe{
      \item{\code{\link{manylm}}}{Fitting Linear Models for Multivariate Abundance Data}
      \item{\code{\link{summary.manylm}}}{summarizie Multivariate Linear Model Fits for Abundance Data}
      \item{\code{\link{anova.manylm}}}{obtain ANOVA for Multivariate Linear Model Fits for Abundance Data}
      \item{\code{\link{plot.manylm}}}{plot diagnostics for a \code{manylm} Object}
   }

\bold{For estimating and displaying Generalized Linear Models:}
    \describe{
      \item{\code{\link{manyglm}}}{fit Generalized Linear Models for Multivariate Abundance Data}
      \item{\code{\link{summary.manyglm}}}{summarize Multivariate Generalized Linear Model Fits for Abundance Data}
      \item{\code{\link{anova.manyglm}}}{obtain Analysis of Deviance for Multivariate Generalized Linear Model Fits for Abundance Data}
      \item{\code{\link{plot.manyglm}}}{plot diagnostics for a \code{manyglm} Object}
    }
Other generic functions like \code{residuals}, \code{predict}, \code{AIC} can be applied to \code{\link{manyglm}} objects.

\bold{For estimating and displaying 'fourth corner models'} with species traits as well as environmental predictors:
    \describe{
      \item{\code{\link{traitglm}}}{predict abundance using a GLM as a function of traits as well as environmental variables}
      \item{\code{\link{anova.traitglm}}}{obtain Analysis of Deviance for a fourth corner model of abundance}
    }
Other generic functions like \code{plot}, \code{residuals}, \code{predict}, \code{AIC} can be applied to \code{\link{traitglm}} objects. Note \code{\link{traitglm}} can work slowly, as it fits a single big model to vectorised data (then wants to resample it when you call \code{anova.traitglm}). 

\bold{For fitting more flexible models:}
    \describe{
      \item{\code{\link{manyany}}}{simultaneously fit univariate models to each response variable from 'any' input function}
      \item{\code{\link{anova.manyany}}}{simultaneously test for a community-level effect, comparing two or more \code{manyany} objects}
      \item{\code{\link{glm1path}}}{fit a path of Generalised Linear Models with L1 ('LASSO') penalties}
      \item{\code{\link{cv.glm1path}}}{choose the value of the L1 penalty in a \code{glm1path} fit by cross-validation}
  }
Other generic functions like \code{residuals}, \code{predict}, \code{AIC} can be applied to \code{\link{manyany}} and \code{\link{glm1path}} objects. These functions also can be on the slow side, especially if all rare species are included.

\bold{For providing a data structure:}
    \describe{
      \item{\code{\link{mvabund}}}{create a mvabund object}
      \item{\code{\link{mvformula}}}{create Model Formulae for Multivariate Abundance Data}
    }
  
\bold{Example datasets:}
    \describe{
      \item{\code{\link{Tasmania}}}{ meiobenthic community data from Tasmania. Used to demonstrate test for interaction. }
      \item{\code{\link{solberg}}}{ solberg species counts with a 3-level treatment factor. }
      \item{\code{\link{spider}}}{ hunting spiders counts from different sites.}
      \item{\code{\link{tikus}}}{ solberg nematode counts from Tikus island.}
      \item{\code{\link{antTraits}}}{ ant counts from Eucalypt forests, with trait measurements.}
    }

For more details, see the documentation for any of the individual functions listed above.

}
\references{

Brown AM, Warton DI, Andrew NR, Binns M, Cassis G and Gibb H (2014) The fourth corner solution - using species traits to better understand how species traits interact with their environment, \emph{Methods in Ecology and Evolution} 5, 344-352.

Warton D.I. (2008a). Raw data graphing: an informative but under-utilized tool for the analysis of multivariate abundances. \emph{Austral Ecology} 33, 290-300.

Warton D.I. (2008b). Penalized normal likelihood and ridge regularization of correlation and covariance matrices. \emph{Journal of the American Statistical Association} 103, 340-349.

Warton D.I. (2011). Regularized sandwich estimators for analysis of high dimensional data using generalized estimating equations. \emph{Biometrics}, 67, 116-123.

Warton DI, Shipley B & Hastie T (2015) CATS regression - a model-based approach to studying trait-based community assembly, \emph{Methods in Ecology and Evolution} 6, 389-398.

Warton D. I., Wright S., and Wang, Y. (2012). Distance-based multivariate analyses confound location and dispersion effects. \emph{Methods in Ecology and Evolution}, 3, 89-101.

Wang Y., Neuman U., Wright S. and Warton D. I. (2012). mvabund: an R package
for model-based analysis of multivariate abundance data. \emph{Methods in Ecology and Evolution}, 3, 471-473. 

}

\author{David Warton \email{David.Warton@unsw.edu.au}, Yi Wang and Ulrike Naumann.}

\seealso{ 
\code{\link{plot.mvabund}}, \code{\link{meanvar.plot}},
\code{\link{manyany}}, \code{\link{manylm}}, \code{\link{manyglm}}, \code{\link{traitglm}}, \code{\link{summary.manylm}}, \code{\link{anova.manyany}}, \code{\link{anova.manylm}}, \code{\link{anova.traitglm}}, \code{\link{anova.manyglm}}, \code{\link{plot.manylm}}
}

\examples{
require(graphics)

## Load the spider dataset:
data(spider)

## Create the mvabund object spiddat:
spiddat <- mvabund(spider$abund)
X <- as.matrix(spider$x)

## Draw a plot of the spider data:
plot(spiddat, col="gray1", n.vars=8, transformation="sqrt", 
xlab=c("Hunting Spider"), ylab="Spider Species", scale.lab="s",
t.lab="t", shift=TRUE, fg= "lightblue", col.main="red", main="Spiders") 


## A mean-variance plot, data organised by year, 
## for 1981 and 1983 only, as in Figure 7a of Warton (2008a):
data(tikus)
tikusdat <- mvabund(tikus$abund)
year <- tikus$x[,1]
is81or83 <- year==81 | year==83
meanvar.plot(tikusdat~year,legend=TRUE, subset=is81or83, col=c(1,10)) 	

## Create a formula for multivariate abundance data:
foo <- mvformula( spiddat~X )

## Create a List of Univariate Formulas:
fooUni <- formulaUnimva(spiddat~X)
fooUniInt <- formulaUnimva(spiddat~X, intercept=TRUE)

## Find the three variables that best explain the response:
best.r.sq( foo, n.xvars= 3)

## Fit a multivariate linear model:
foo <- mvformula( spiddat~X )
lm.spider <- manylm(foo)

## Plot Diagnostics for a multivariate linear model:
plot(lm.spider,which=1:2,col.main="red",cex=3,overlay=FALSE)

## Obtain a summary of test statistics using residual resampling:
summary(lm.spider, nBoot=500)

## Calculate a ANOVA Table:
anova(lm.spider, nBoot=500)

}
    
\docType{ package }
\keyword{ documentation }
