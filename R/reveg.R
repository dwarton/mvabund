#' Invertebrate abundances in a revegetation study
#'
#' Data from a study looking at the effect of revegetation on invertebrate communities 
#' (data from Anthony Pik, Macquarie University). Invertebrates were sampled in 4-5 pitfall
#' traps at eight sites that had undergone revegetation, and two sites that hadn't, and it
#' was of interest to see what the effect of revegetation was on the invertebrate community.
#'
#' @docType data
#'
#' @usage data(reveg)
#'
#' @format A list containing three objects:\describe{
#' \item{abund}{A data frame with abundances of 24 Orders of invertebrate.}
#' \item{pitfalls}{A vector specifying the number of pitfall traps that were collected at each site.}
#' \item{treatment}{Whether the site was a `Control` or a site that had undergone revegetation (`Impact`).}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(reveg)
#' worms = reveg$abund$Haplotaxida
#' plot(worms~treatment, data=reveg, horizontal=TRUE,
#'   las=1, xlab="",ylab="Worm abundance")
#' 
"reveg"