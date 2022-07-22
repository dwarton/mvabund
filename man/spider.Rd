\name{spider}
\docType{data}
\alias{spider}
\title{Spider data}
\description{
data from spider2 directory, CANOCO FORTRAN package, with trait variables added.
}
\usage{
data(spider)
}
\format{
A list containing the elements
\describe{
\item{abund}{
A matrix with 28 observations of abundance of 12 hunting spider species}
\item{x}{
A data frame of six (transformed) environmental variables at each of the 28 sites.
}
\item{trait}{
A data frame of three species traits for each of the 12 species.
}
}
The matrix \code{abund} has the following variables
   \describe{
    \item{Alopacce}{(numeric) Abundance of the species Alopecosa accentuata }
    \item{Alopcune}{(numeric) Abundance of the species Alopecosa cuneata }
    \item{Alopfabr}{(numeric) Abundance of the species Alopecosa fabrilis }
    \item{Arctlute}{(numeric) Abundance of the species Arctosa lutetiana }
    \item{Arctperi}{(numeric) Abundance of the species Arctosa perita }
    \item{Auloalbi}{(numeric) Abundance of the species Aulonia albimana }
    \item{Pardlugu}{(numeric) Abundance of the species Pardosa lugubris }
    \item{Pardmont}{(numeric) Abundance of the species Pardosa monticola }
    \item{Pardnigr}{(numeric) Abundance of the species Pardosa nigriceps }
    \item{Pardpull}{(numeric) Abundance of the species Pardosa pullata }
    \item{Trocterr}{(numeric) Abundance of the species Trochosa terricola }
    \item{Zoraspin}{(numeric) Abundance of the species Zora spinimana }
   }

The data frame \code{x} has the following variables
   \describe{
    \item{soil.dry}{(numeric) Soil dry mass }
    \item{bare.sand}{(numeric) Cover bare sand }
    \item{fallen.leaves}{(numeric) Cover fallen leaves / twigs }
    \item{moss}{(numeric) Cover moss}
    \item{herb.layer}{(numeric) Cover herb layer}
    \item{reflection}{(numeric) Reflection of the soil surface with a cloudless sky}
   }
   These variables have already been log(x+1)-transformed.

The data frame \code{trait} was constructed by Googling each species and recording variables from species descriptions and images of specimens:
   \describe{
    \item{length}{(numeric) Length (log-transformed), averaged across typical lengths (in centimetres) for male and females}
    \item{colour}{(factor) Predominant colour, "yellow" or "dark"}
    \item{marks}{(factor) Whether the spider typically has markings on it: "none", "spots" or "stripes"}
   }
}

\source{Data attributed to van der Aart & Smeenk-Enserink (1975), obtained from the spider2 directory, CANOCO FORTRAN package. Trait data largely extracted from Wikipedia entries for the species.}
\references{

ter Braak, C. J. F. and Smilauer, P.  (1998)  CANOCO reference manual and user's guide to CANOCO for Windows: software for canonical community ordination (version 4). Microcomputer Power, New York, New York, USA.

van der Aart, P. J. M., and Smeenk-Enserink, N. (1975) Correlations between
distributions of hunting spiders (Lycos- idae, Ctenidae) and environmental
characteristics in a dune area. \emph{Netherlands Journal of Zoology} \bold{25},
1-45.
}
\details{The abundance of each species was measured as a count of the number of organisms in the sample.}
%\seealso{}
\examples{
require(graphics)

data(spider)
spiddat <- as.mvabund(spider$abund)

plot(spiddat)
}
\keyword{datasets}
