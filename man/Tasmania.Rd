\name{Tasmania}
\docType{data}
\alias{Tasmania}
\title{Tasmania Dataset}
\description{
This dataset contains a list with community abundance data of species and two factor variables, namely treatment and block. See (Warwick et.al. (1990)) for more details.
}
\usage{
data(Tasmania)
}
\format{
A list containing the elements
\describe{
        \item{abund}{A data frame with 16 observations of 56 Meiobenthos species exposed to a disturbance treatment in a spatially blocked design. Four blocks of four samples were collected such that each block comprised of two disturbed and undisturbed samples.}
        \item{copepods}{A subset of \code{abund} of 12 Copepod species.}
        \item{nematodes}{A subset of \code{abund} of 39 Nematode species.}
        \item{treatment}{A two-level factor veraible.}
        \item{block}{A four-level factor variable.}
}
}

\details{
The count data (number of each Meiobenthos species in each sample) were collected in a spatia
lly blocked design. The labels are made to the four replicate cores within each block,
with \code{B} labeling for the block ID and \code{D} labeling for the disturbed sample ID and \code{U} labeling for the undisturbed sample ID. The data frame \code{abund} contains 12 Copepod species, 39 Nematode species and 4 undetermined ones.

The 12 Copepod species are: \cr
Ameira,      Adopsyllus,  Ectinosoma,  Ectinosomat, Haloschizo,\cr
Lepta.A,     Lepta.B,     Lepta.C,     Mictyricola, Parevansula,\cr
Quin,        Rhizothrix\cr

The 39 Nematode species are: \cr
Actinonema,        Axonolaimus,       Bathylaimus, \cr
Calyptronema,      Chaetonema,        Chromaspirina, \cr
Comesoma,          Daptonema,         Desmodora.A, \cr
Desmodora.B,       Enoploides,        Enoplus, \cr
Epacanthion.A,     Epacanthion.B,     Eubostrichus, \cr
Eurystomina,       Hypodontolaimus.A, Hypodontolaimus.B, \cr
Leptolaimus,       Leptonemella,      Mesacanthion, \cr
Microlaimus,       Monhystera,        Nannoluimoides.A, \cr
Nannolaimoides.B,  Neochromadora.A,   Neochromadora.B, \cr
Odontophora,       Oncholaimus,       Qnvx, \cr
Paracanthonchus,   Polysigma,         Praeacanthenchus, \cr
Promonhystera,     Pseudosteineria,   Sabatieria, \cr
Spilophorella,     Symplocostoma,     Viscosia \cr

The data frame \code{copepod} stores the subset of 12 Copepod species, and the data frame \code{nematode} stores the subset of 39 Nematode species.

\code{treatment} indicates disturbed or undisturbed treatment for the 16 observations of each species in the Tasmania dataset.

\code{block} indicates the block ID for the 16 observations of each species in the Tasmania dataset.
}

\references{
Warwick, R. M., Clarke, K. R. and Gee, J. M. (1990). The effect of disturbance by soldier crabs Mictyris platycheles H. Milne Edwards on meiobenthic communiy structure. \emph{J. Exp. Mar. Biol. Ecol.}, \bold{135}, 19-33.
}
\examples{
require(graphics)

data(Tasmania)
tasm.cop <- mvabund(Tasmania$copepods)
treatment <- Tasmania$treatment
block <- Tasmania$block

plot(tasm.cop~block*treatment)

}
\keyword{datasets}
