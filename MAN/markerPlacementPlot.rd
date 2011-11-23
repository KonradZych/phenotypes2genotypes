\name{markerPlacementPlot}
\alias{markerPlacementPlot}


\title{Plot number of markers selected.}

\description{
 Plot number of markers selected with different thresholds.
}

\usage{
	markerPlacementPlot(population, cross)
}

\arguments{
 \item{population}{ an object of class \code{\link{population}}}
 \item{cross}{ an object of R/qtl class cross}
}

\value{
  Plot.
}

\details{
This plot is really usefull while saturating existing map. It helps choose best threshold for marker selection, showing how much markers will
be selected with different threshold values.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastCross)
	data(yeastPopulation)
	assignment <- createNewMap(yeastPopulation,yeastCross,n.chr=16,verbose=TRUE,map="physical",comparisonMethod=sumMajorityCorrelation, use.orderMarkers=FALSE,reOrder=FALSE)
  assignment #boring,but expected
  ordering <- assignedChrToMarkers(assignment,yeastCross)
}

\seealso{
  \code{\link{saturateExistingMap}} - saturate existing map
}

\keyword{manip}
