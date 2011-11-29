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
	markerPlacementPlot(yeastPopulation,yeastCross)
}

\seealso{
  \code{\link{cross.saturate}} - saturate existing map
}

\keyword{manip}
