\name{plot.mvabund}
\alias{plot.mvabund}
\alias{plot.mvformula}
\title{Plot Multivariate Abundance Data and Formulae}
\description{
Produces a range of plots for visualising multivariate abundance data and 
its relationship to environmental variables, including: dot-plots and boxplots 
for different levels of a factor
stacked by response variable; comparative dot-plots and boxplots for different 
levels of a factor, stacked by response variable; scatterplots of abundances 
against a set of explanatory variables; scatterplots of pair-wise abundance 
data, e.g. from repeated measures. See details below.
}
\usage{
\method{plot}{mvabund}(x, y, type="p", overall.main="", n.vars=12, 
  var.subset=NA, transformation="log", ...)

\method{plot}{mvformula}(x,y=NA, type="p", var.subset=NA, 
  n.vars= if(any(is.na(list(var.subset)))) 12 else length(var.subset),
  xvar.select=TRUE, xvar.subset = NA, n.xvars=NA, transformation="log", ...) 
}


\arguments{
\item{x}{
for the \code{mvabund} method, x is a \code{mvabund} object.\cr
For the \code{mvformula} method, x is a \code{mvformula} object,
a Model Formula to be used.
}
\item{y}{in \code{plot.mvabund} an optional second matrix with multivariate
abundance data in \code{plot.mvformula} an optional matrix of the
independent variables to explain x.}
\item{type}{what type of plot should be drawn. Useful types are "p" for
scatterplot, "bx" for boxplot and "n" for no plot. Other types, see
\code{plot} are allowed, but usually don't give a meaningful output.}
\item{overall.main}{a character to display as title for every window.}
\item{var.subset}{a numeric vector of indices indicating which variables of the
mvabund.object should be included on the plot.}
\item{n.vars}{the number of variables to include in the plot.}
\item{xvar.select}{whether only a subset of x variables should be plotted or all.}
\item{n.xvars}{the number of the most relevant x variables that should be plotted,
 is not used if \code{xvar.select = FALSE}. If NA it will be set to at most 3.}
\item{xvar.subset}{a subset of x variables that should be plotted, is not used if \code{xvar.select = FALSE}.}
\item{transformation}{an optional transformation, if no formula is given,
"no" = untransformed, "sqrt"=square root transformed,
"log" (default)=log(Y/min+1) transformed, "sqrt4" =4th root transformed.\cr
Note that if \code{plot.mvabund} is called explicitly, and two data objects 
supplied, none of which is a \code{mvabund} object, then \code{plot.mvformula}
will be called (See Details). The argument \code{transformation} is then NOT 
available. }
\item{\dots}{arguments to be passed to or from other methods.}
}
\details{
The function \code{plot.mvabund} produces plots for the visualisation of
multivariate abundance data and their relationships to environmental variables.
The approach taken is to separately plot the relationship between each response
variable and environmental
variables, that is, to visualise the marginal distribution, as in Warton (2008). 
Three main types of plot that can be produced:

(1) Dot-plots or boxplots stacked along the y-axis by response variable. If a
factor is given, comparative dot-plots/boxplots are produced, comparing
abundances across each factor level. This type of plot is produced when one
multivariate abundance dataset is given as an input argument, either on its own,
or together with a factor, as in the examples using the solberg dataset below.

(2) Scatterplots of multivariate abundances against environmental variables,
with separate plots for separate response variables. This type of plot is
produced when one multivariate abundance dataset is given as an input argument
together with an environmental variable or a set of environmental variables.

(3) Scatterplots of a paired sample of multivariate abundances. This type of
plot is produced when two multivariate abundance datasets are given as input
arguments, and their size and variable names match each other. It is up to the
user to ensure that the rows match for these two datasets.

There are several methods for calling \code{plot.mvabund}:  \cr
(a) plot.mvabund("matrix", ...) \cr
The multivariate abundances are stored
in the data matrix. Including an optional second argument determines
whether a plot of type (1) is produced (if no second argument or if it is a factor), 
or a plot of type (2) (if one or a set of environmental variables is given), or a 
plot of type (3) (if a second matching multivariate abundance dataset is given). \cr
Instead of a matrix, \code{mvabund} objects can be used. \cr
(b) plot("mvabund object", ...) \cr
You can define mvabund objects using the function 
\code{mvabund}. Then the behaviour of the plot function is the same as 
\code{plot.mvabund} above. \cr
(c) plot.mvformula("response"~"terms") \cr
The first of these two objects must be the multivariate abundances, which can be
either a matrix or a \code{mvabund} object. The terms determine the type of
plot produced. The terms can be either a single vector or matrix or 
a number of vectors or matrices, separated by \code{+}. 
Compare \code{formula} for further details on the specification
of the terms. \cr
(d) plot("mvformula object") \cr
You can define mvformula objects using the function \code{mvformula}.
Note that the response cannot be a data frame object.

For plots of type (3) above, you must use method (a) or (b).
Plot methods (c) and (d) require that both the response and explanatory
variables are specified, i.e. formulas like '~x' or 'y~1' cannot be plotted.

See below examples to see how each of these methods is applied.

Multivariate abundance datasets typically have many variables,
more than can be visualised in a single window, so by default plot.mvabund
subsets abundance variables (and where appropriate, environmental variables).
By default the 12 most abundant variables are plotted (determined on transformed
 variables if the response is transformed in the mvformula method),
although this setting can be controlled via the argument \code{n.vars}, and the
variables included in the subset to be plotted can be controlled via
\code{var.subset}. It is possible for example to plot the abundance variables
most significantly associated with an environmental variable,
as in the Solberg example below.

To produce boxplots rather than dot-plots in type (1) plots, use the argument
\code{type="bx"}.

For type (2) plots, if only one environmental variable is specified, plots for 
different abundance variables are arranged in a rectangular array (different 
abundance variables in different rows and columns). If however more than one
environmental variable is specified, different columns correspond to different
environmental variables (and different abundance variables in different rows).
If more than 3 environmental variables are specified, the 3 will be selected
that maximise average R^2 when \code{manylm} is applied (using the subset
selection function \code{best.r.sq}). To switch off this subset selection, set
\code{xvar.select=FALSE}, or choose your own subset of environmental variables
using \code{xvar.subset}.

To control the appearance of points on dot-plots and scatterplots, usual
arguments apply (see \code{par} for details). The plotting symbols \code{pch}
and their \code{color} can be a vector, and if the plot function is called via
a mvformula object, it can also be a list, where the list elements corresponds
to the symbols / colors used in the plots for different
independent variables. 

If some of the formula terms are factor variables, these will be drawn in
boxplots.
Note, that the plots produced by \code{plot.mvformula}
depend on whether the first independent variable is a factor or not. 
See the examples for the different possibilities of boxplots that can be
produced.

If two objects are passed and only one of them is an \code{mvabund} object,
the resulting plots will be the same as if a formula was supplied, having the
\code{mvabund} object as response variable.\cr
If both objects are not \code{mvabund} objects, it will be tried to guess which
one of them is the response. The following logic applies:
If \code{y} is not a \code{data.frame}, it will be assumed that \code{y} is the
response. Note that \code{y} is the second object, if argument names are not
supplied.
If \code{y} is a \code{data.frame} and \code{x} is not a \code{data.frame},
it will be assumed that \code{x} is the response. If both objects are
data frames an error will result, as the function is designed for \code{mvabund}
objects!

The argument \code{shift} controls whether or not points are shifted on dotplots 
so that they do not overlap. This argument is ignored for boxplots and
scatterplots (type (2) or type (3) graphs).
}
\section{Warning}{
The argument \code{log}, that is available in lots of plotting functions can not
be used for plotting \code{mvabund} or \code{mvformula} objects. Instead use
\code{transformation} for the \code{mvabund} method and for the
\code{mvformula} method include any transformations in the formula.
}
\references{
Warton D.I. (2008). Raw data graphing: an informative but under-utilized tool
for the analysis of multivariate abundances. \emph{Austral Ecology} 33(3), 290-300.
}
\author{
Ulrike Naumann, Yi Wang, Stephen Wright and David Warton <David.Warton@unsw.edu.au>.
}
\seealso{
\code{\link{boxplot.mvabund}}, \code{\link{meanvar.plot}}, \code{\link{plot.manylm}}, 
\code{\link{plot.manyglm}}.
}
\examples{

require(graphics)

############################
## Some "type (1)" plots ##
############################

data(solberg)
solbdat <- solberg$abund
treatment<- solberg$x

## Plot a multivariate dataset (Species vs Abundance)
plot.mvabund(solbdat)

## Alternatively, the plot command could be used directly if spiddat is 
## defined as an mvabund object:
solbmvabund <- mvabund(solbdat)
plot(solbmvabund)

## Draw an mvabund object in a boxplot, but using the 20 most abundant 
## variables in the plot, using the square root transform, and adding 
## coloured axes and title:
plot.mvabund(solbdat, n.vars=20, type="bx", transformation="sqrt", 
fg="lightblue", main="Solberg abundances", col.main="red")

## Plot Species (split by treatment) vs Abundance
plot(solbmvabund,treatment)

## This can also be produced using
plot(solbmvabund~treatment)

## To use plot.mvabund to plot only the variables with P-values less than 0.1:
lm.solberg <- manylm(log(solbmvabund+1)~treatment)
anova.solb <- anova(lm.solberg, p.uni="unadjusted")
pj = anova.solb$uni.p

plot(solbmvabund~treatment, var.subset=pj<0.1)

## Or to plot only the 12 most significant variables, according to their 
## univariate ANOVA P-values:
pj.sort = sort(pj, index.return=TRUE)
plot(solbmvabund~treatment, var.subset=pj.sort$ix[1:12])


############################
## Some "type (2)" plots ##
############################
#load and convert data
data(spider)
spiddat <- mvabund(spider$abund)
spidx <- mvabund(spider$x)

#create labels vectors
pch.vec <- as.numeric(spidx[,3]<2)
pch.vec[pch.vec==0] <- 3

#Scale the soil water variable
soilWater <- spidx[,1]

#Create the Table for the main titles of each plot
title <- c("\n\nAlopecosa accentuata", "\n\nAlopecosa cuneata",
           "\n\nAlopecosa fabrilis", "\n\nArctosa lutetiana", 
           "\n\nArctosa perita", "\n\nAulonia albimana", 
           "\n\nPardosa lugubris", "\n\nPardosa monticola", 
           "\n\nPardosa nigriceps", "\n\nPardosa pullata",
           "\n\nTrochosa terricola", "\n\nZora spinimana")

#Plot Species Abundance vs Environmental variable
plot.mvformula(log(spiddat+1) ~ exp(soilWater), main=title, 
xlab="\% Soil Moist - Log Scale  ", ylab="Abundance [log scale]", 
overall.main="Species Abundance vs \%Soil Moisture", col=pch.vec, 
fg="grey", pch=pch.vec, las=1, scale.lab="ss",t.lab="o", mfrow=c(4,3),log="x")

#Add a Margin
par(xpd=NA)
legend("topright",pch=c(1,3),col=c(1,3),legend = c("few twigs", "many twigs"), 
cex = 1, inset=c(0,-0.19))


############################
## Some "type (3)" plots ##
############################

##Plot 1981 Abundance vs 1983 Abundance
data(tikus)
year <- tikus$x[,1]
tikusdat <- mvabund(tikus$abund)
site <- tikus$x[,2]

plot(tikusdat[year==81,], tikusdat[year==83,], col.main="blue", 
xlab="1981 abundance", ylab="1983 abundance")

}
\keyword{hplot}

