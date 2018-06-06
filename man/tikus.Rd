\name{tikus}
\docType{data}
\alias{tikus}
\title{Tikus Island Dataset}
\description{
This dataset contains a list with abundance data of species at different locations
in the Tikus island and explanatory variables.
}
\usage{
data(tikus)
}
\format{
A list containing the elements
\describe{
\item{abund}{A data frame with 60 observations at different locations
 of abundances on 75 variables, the species. See Details.
}

\item{x}{A data frame containing the id information for the Tikus island dataset.
The data frame has 60 observations on 2 variables. See Details.
}
}
}
\details{
The abundance of each species was measured as the length (in centimetres) of a
10 metre transect which intersected with the species.

\code{tikus} is a list containing the elements \code{abund} and \code{x}.
The data frame \code{abund} contains 75 variables, the species:

Psammocora contigua, Psammocora digitata, Pocillopora damicornis,
Pocillopora verrucosa, Stylopora pistillata, Acropora bruegemanni,
Acropora robusta, Acropora grandis, Acropora intermedia,
Acropora formosa, Acropora splendida, Acropera aspera,
Acropora hyacinthus, Acropora palifera, Acropora cytherea,
Acropora tenuis, Acropora pulchra, Acropora nasuta,
Acropora humilis, Acropora diversa, Acropora digitifera,
Acropora divaricata, Acropora subglabra, Acropora cerealis,
Acropora valida, Acropora acuminata, Acropora elsevi,
Acropora millepora, Montipora monasteriata, Montipora tuberculosa,
Montipora hispida, Montipora digitata, Montipora foliosa,
Montipora verrucosa, Fungia fungites, Fungia paumotensis,
Fungia concina, Fungia scutaria, Halomitra limax, Pavona varians,
Pavona venosa, Pavona cactus, Coeloseris mayeri,
Galaxea fascicularis, Symphyllia radians, Lobophyllia corymbosa,
Lobophyllia hemprichii, Porites cylindrica, Porites lichen,
Porites lobata, Porites lutea, Porites nigrescens, Porites solida,
Porites stephensoni, Goniopora lobata, Favia pallida, Favia speciosa,
Favia stelligera, Favia rotumana, Favites abdita, Favites chinensis,
Goniastrea rectiformis, Goniastrea pectinata, Goniastrea sp,
Dulophyllia crispa, Platygyra daedalea, Platygyra sinensis,
Hydnopora rigida, Leptastrea purpurea, Leptastrea pruinosa,
Cyphastrea serailia, Millepora platyphylla, Millepora dichotoma,
Millepora intrincata, Heliopora coerulea

\code{x} has the following variables:
   \describe{
    \item{time}{(factor) the year in which the measurement was taken.}
    \item{rep}{(factor) the location id. }
   }
}

\references{
R.M. Warwick, K.R. Clarke and Suharsono (1990) A statistical analysis of coral
community responses to the 19823 El Nino in the Thousand Islands, Indonesia,
\emph{Coral Reefs} \bold{8}, 171179.
}
\examples{
require(graphics)

data(tikus)

tikusdat <- as.mvabund(tikus$abund)
tikusid	<- tikus$x
foo	<- mvformula(tikusdat~tikusid[,1] + tikusid[,2])

plot(foo)
}
\keyword{datasets}
