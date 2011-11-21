\name{markersCorPlot}
\alias{markersCorPlot}

\title{Plotting correlation between markers on two maps.}

\description{
  Plotting correlation between two maps together with markers placement (comparison of coverage).
}

\usage{
	markersCorPlot(cross, population, map=c("genetic","physical"), cmBetween=25, comparisonMethod = c(sumMajorityCorrelation,majorityCorrelation,meanCorrelation), chr, verbose=TRUE)
}

\arguments{
 \item{cross}{ R/qtl cross type object.}
  \item{population}{ an object of class \code{\link{population}}}
 \item{map}{ which map (from ones stored in population$maps) should be used fo assigning chromosomes on the created map}
 \item{cmBetween}{ offset between chromosomes specified in cM}
 \item{comparisonMethod}{method used tocompare chromosomes from the new map to the original ones while assigning
   \itemize{
    \item{\code{\link{sumMajorityCorrelation}}}{}
    \item{\code{\link{majorityCorrelation}}}{}
    \item{\code{\link{meanCorrelation}}}{}
    \item{\code{\link{majorityOfMarkers}}}{}
  }
 }
 \item{chr}{ specifies subset of chromosomes to be shown }
 \item{verbose}{ be verbose}
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
	markersCorPlot(yeastCross,yeastPopulation,map="physical")
}

\seealso{
  \code{\link{plotMapComparison}} - Plotting routine for comparison of two genetic maps.
  \code{\link{projectOldMarkers}} - Plotting routine for showing how markers from original map are placed on saturated map.
}

\keyword{manip}
