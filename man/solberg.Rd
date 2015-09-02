\name{solberg}
\docType{data}
\alias{solberg}
% \encoding{latin1}
\title{Solberg Data}
\description{
This dataset contains a list with abundance data of species and a factor variable.
}
\usage{
data(solberg)
}
\format{
A list containing the elements

\describe{
\item{abund}{a data frame containing 12 rows and has 53 variables, corresponding to the species. It has the following variables:\cr
\code{Paramesacanthion_sp.}, \code{Halalaimus_sp.}, \code{Viscosia_sp.}, \cr
\code{Symplocostoma_sp.}, \code{Bathylaimus_inermis}, \code{Bathylaimus_sp.},\cr
\code{Rhabdodemania_sp.}, \code{Pandolaimus_latilaimus}, \code{Halanonchus_sp.}
,\cr
\code{Trefusia_sp.}, \code{Chromadora_sp.}, \code{Dichromadora_sp.},\cr
\code{Neochromadora_sp.}, \code{Prochromadorella_sp.}, \code{Neotonchus_sp.},\cr
\code{Marylynnia_complexa}, \code{Paracanthonchus_sp.}, \code{Cyatholaimidae_un
.},\cr
\code{Choniolaimus_papillatus}, \code{Halichoanolaimus_dolichurus},\cr
\code{Richtersia_inaequalis}, \code{Dorylaimopsis_punctatus},\cr
\code{Sabatieria_longicaudata}, \code{Sabatieria_punctata},\cr
\code{Sabatieria_sp.}, \code{Setosabieria_hilarula},\cr
\code{Chromaspirina_sp.}, \code{Molgolaimus_sp.}, \cr
\code{Spirinia_parasitifera}, \code{Aponema_torosa}, \cr
\code{Microlaimus_sp.1}, \code{Microlaimus_sp.2}, \code{Camacolaimus_sp.}, \cr
\code{Leptolaimus_elegans}, \code{Monhystera_sp.},\cr 
\code{Amphimonhystera_sp.}, \code{Daptonema_sp.1}, \code{Daptonema_sp.2}, \cr
\code{Daptonema_sp.3}, \code{Theristus_aff_acer}, \code{Xyalidae_un.},\cr
\code{Sphaerolaimus_macrocirculus}, \code{Sphaerolaimus_paradoxus},\cr
\code{Desmolaimus_sp.}, \code{Eleutherolaimus_sp.}, \code{Eumorpholaimus_sp.},\cr
\code{Terschellingia_longicaudata}, \code{Paralinhomoeus_conicaudatus},\cr
\code{Linhomieidae_un.A}, \code{Linhomieidae_un.B}, \code{Axonolaimus_sp.},\cr
\code{Odontophora_sp.}, \code{Unidentified}
}
\item{x}{a factor with the levels \code{control}, \code{high}, \code{low}}
}
}
%\source{}
\references{
Gee J. M., Warwick R. M., Schaanning M., Berge J. A. and Ambrose W. G. Jr (1985) Effects of organic enrichment on meiofaunal abundance and community structure in sublittoral soft sediments. \emph{Journal of Experimental Marine Biology and Ecology}. 91(3), 247-262.

}
\details{The abundance of each species was measured as the count 
of the number of organisms in the sample.}
% \seealso{}
\examples{
data(solberg)
solbergdat <- mvabund( solberg$abund )
treatment <- solberg$x

## Create a formula for multivariate abundance data:
foo.sol <- mvformula( solbergdat ~ treatment )

## Fit a multivariate linear model:
lm.solberg <- manylm(foo.sol)
lm.solberg
}
\keyword{datasets}
