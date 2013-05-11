\name{qtl.comparison.plot}
\alias{qtl.comparison.plot}

\title{Comparison of qtl profiles.}

\description{
  Plots comparison between the qtl profiles of two cross objects.
}

\usage{
	qtl.comparison.plot(cross1, cross2, chr, ...)
}

\arguments{
 \item{cross1}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details.}
 \item{cross2}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{ Specifies the chromosome to be shown (only one chromosome can be plotted at a time.}
 \item{...}{ Arguments passed to scanone function (see \code{\link[qtl]{scanone}}).}
}

\details{
Plots markers from moth old and new map as points and in the background - comparison between them done using selected comparison method.
}

\value{
	Matrix of comparisons between chromosomes obtained using comparison method.
}

\author{
	Konrad Zych \email{k.zych@rug.nl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{k.zych@rug.nl}
}

\examples{
	data(yeastCross)
	qtl.comparison.plot(yeastCross,yeastCross)
}

\seealso{
  \itemize{
    \item{\code{\link{plotMapComparison}}}{ -  Plotting routine for comparison of two genetic maps.}
    \item{\code{\link{projectOldMarkers}}}{ -   Plotting routine for showing how markers from original map are placed on saturated map.}
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or chromosome assignment vector.}
}
}

\keyword{manip}
