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
 \item{cross}{ R/qtl cross type object.}
  \item{population}{ an object of class \code{\link{population}}}
 \item{map}{ which map (from ones stored in population$maps) should be used fo assigning chromosomes on the created map}
 \item{chr}{ specifies subset of chromosomes to be shown }
}

\value{
	Plot.
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
  \code{\link{projectOldMarkers}} - Plotting routine for showing how markers from original map are placed on saturated map.
  \code{\link{markersCorPlot}} - Plotting correlation between two maps together with markers placement (comparison of coverage).
}

\keyword{manip}
