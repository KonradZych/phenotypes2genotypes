\name{projectOldMarkers}
\alias{projectOldMarkers}

\title{Projecting old markers on saturated map.}

\description{
  Plotting routine for showing how markers from original map are placed on saturated map.
}

\usage{
	projectOldMarkers(cross,population,map=c("genetic","physical"),label=c("positions","names"))
}

\arguments{
 \item{cross}{ R/qtl cross type object.}
  \item{population}{ an object of class \code{\link{population}}}
 \item{map}{ which map (from ones stored in population$maps) should be used fo assigning chromosomes on the created map}
 \item{label}{ how old markers should be labeled on the plot, with their positions or their names}
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
	projectOldMarkers(yeastCross,yeastPopulation,map="physical")
}

\seealso{
  \code{\link{plotMapComparison}} - Plotting routine for comparison of two genetic maps.
  \code{\link{markersCorPlot}} - Plotting correlation between two maps together with markers placement (comparison of coverage).
}

\keyword{manip}
