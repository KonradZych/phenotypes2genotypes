\name{comparisonMethod}
\alias{comparisonMethod}
\alias{sumMajorityCorrelation}
\alias{majorityCorrelation}
\alias{meanCorrelation}
\alias{majorityOfMarkers}

\title{Methods for comparing maps.}

\description{
  Different methods for comparing maps
}

\usage{
	sumMajorityCorrelation(cross, originalMap, population, verbose=FALSE)
	majorityCorrelation(cross, originalMap, population, verbose=FALSE)
	meanCorrelation(cross, originalMap, population, verbose=FALSE)
	majorityOfMarkers(cross, originalMap, population, verbose=FALSE)
}

\arguments{
 \item{cross}{ an object of class cross}
 \item{originalMap}{ an object containing original map used for comparison}
 \item{population}{ an object of class population}
 \item{verbose}{ be verbose}
}

\details{
Using cross.denovo or markersCorPlot functions, user will need choose one of these methods of comparing created map to the original one.
They are described below. However, those functions are always called inside the function, so only description without examples of usage,
is given. 

sumMajorityCorrelation - for every chromosome in cross for every marker checks the marker it is
   having highest correlation with. Checks on which chromosome this marker is placed in old map. For each of
   new chromosomes one or more of chromosomes from old map will be represented. Function sums correlations for
   each pair of those and for every new chromosomes assigns old chromosome with highest cumulative cor.

majorityCorrelation - for every chromosome in cross for every marker checks the marker it is
   having highest correlation with. Checks on which chromosome this marker is placed in old map. For each of
   new chromosomee, old chromosome with most markers with high correlation is assigned:

meanCorrelation - assigning chromosome from new map to old ones using sum of the mean correlation between their markers.

majorityOfMarkers - for each chromosome in the cross object (either created inside the function or provided
  by user) chromosome from original map, where most markers from new chromosome are is assigned.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\seealso{
  \code{\link{reorganizeMarkersWithin}} - apply new ordering on the cross object usign ordering vector
  \code{\link{assignedChrToMarkers}} - create ordering vector from chromosome assignment vector
  \code{\link{cross.denovo}} - creating de novo genetic map or chromosome assignment vector
  \code{\link{reduceChromosomesNumber}} - number of routines to reduce number of chromosomes of cross object
  \code{\link{markersCorPlot}} - Plotting correlation between two maps together with markers placement (comparison of coverage).
}

\keyword{manip}
