\name{plotMapComparison}
\alias{plotMapComparison}
\alias{coloringMode}

\title{Plotting routine for comparison of two genetic maps.}

\description{
  Plotting routine for comparison of two genetic maps.
}

\usage{
	plotMapComparison(cross, population, map=c("genetic","physical"), chr)
}

\arguments{
\item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
\item{population}{ An object of class \code{\link{population}}. See \code{\link{createPopulation}} for details. }
 \item{map}{ Which map (from ones stored in population$maps) should be used fo assigning chromosomes on the created map.}
 \item{chr}{ Specifies subset of chromosomes to be shown.}
}

\value{
	None.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
	data(yeastCross)
	plotMapComparison(yeastCross,yeastPopulation,map="physical")
}

\seealso{
  \itemize{
    \item{\code{\link{markersCorPlot}}}{ -   Plotting correlation between two maps together with markers placement (comparison of coverage).}
    \item{\code{\link{projectOldMarkers}}}{ -   Plotting routine for showing how markers from original map are placed on saturated map.}
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or chromosome assignment vector.}
}
}

\keyword{manip}
