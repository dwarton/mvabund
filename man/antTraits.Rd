\name{antTraits}
\docType{data}
\alias{antTraits}
\title{Ant data, with species traits}
\description{
Abundances of 41 epigaeic ant species across 30 sites in south-eastern Australia, with species trait and environmental data 
}
\usage{
data(antTraits)
}
\format{
A list containing three elements:
\describe{
\item{abund}{A data frame with observations at 30 different locations
 of abundances of 41 epigaeic ant species.
}
\item{env}{A data frame containing 7 environmental variables from transects at each of the 30 sites:
\describe{
\item{Bare.ground}{Percent cover of bare ground, as estimated from ten 1x1 metre quadrats}
\item{Canopy.cover}{Percent canopy cover, as estimated from two 20x20m transects}
\item{Shrub.cover}{Percent canopy cover, as estimated from two 20x20m transects}
\item{Volume.lying.CWD}{Estimated volume of Coarse Woody Debris in two 20x20m transects, including all debris >5cm diameter.}
\item{Feral.mammal.dung}{Proportion of quadrats including mammal dung, out of ten 1x1m quadrats.}
}
}
\item{traits}{A data frame containing 5 species traits measured for each of the 41 species. Weber's length was log-transformed, Femur length was log-transformed then regressed against log(Weber's length), to remove the effect of size.
\describe{
\item{Femur.length}{Residuals from regression of log(Femur length) against log(Weber's length)}
\item{No.spines}{Number of spines on propodeum and petioles, as an integer value}
\item{Pilosity}{A factor with four levels of pilosity, 0 = No or very few hairs; 1 = a sparse but regular covering of hairs; 2 = a consistent, moderate covering of hair; 3 = very dense hair covering}
\item{Polymorphism}{0 = Monomorphic, 1 = polymorphic, 2 = dimorphic}
\item{Webers.length}{log transformed. Body length, measured as the distance from the anterodorsal margin of the pronotum to the posteroventral margin of the propodeum}
}
}
}
}
\references{
Gibb H, Stoklosa J, Warton, DI, Brown, AM, Andrew, NR and Cunningham, SA (2015) Does morphology predict trophic position and habitat use of ant species and assemblages? Oecologia 177, 519-531.
}
\examples{
data(antTraits)
ft = traitglm(antTraits$abund,antTraits$env,antTraits$traits) #to do a fourth corner analysis
}
\keyword{datasets}
