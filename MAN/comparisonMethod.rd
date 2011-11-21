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
	majorityCorrelation(cross, originalMap, population, verbose=FALSE)
	
}

\arguments{
 \item{population}{ an object of class population}
 \item{cross}{ an object of class cross, if not supported, it will be created using data stored in population}
  \item{map}{ 
  Which map should be used for comparison:
  \itemize{
    \item{genetic}{ - genetic map from cross$maps$genetic}
    \item{physical}{ - physical map from cross$maps$physical}
  }
  }
 \item{corTreshold}{ markers not having corelation above this number with any of chromosomes are removed}
 \item{addMarkers}{ should markers from map used for ordering be added to resulting map}
 \item{verbose}{ Be verbose}
}

\value{
  an object of class cross
}

\details{
This function saturates existing map (stored in the population object) with markers derived form gene expression data
(provided inside cross or population. Correlation matrix between those two is made and based on it, new markers sre
being placed on the map.
}
\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
	cross <- saturateExistingMap(yeastPopulation,map="physical",verbose=TRUE,debugMode=2)
}

\seealso{
  \code{\link{reorganizeMarkersWithin}} - apply new ordering on the cross object usign ordering vector
  \code{\link{assignedChrToMarkers}} - create ordering vector from chromosome assignment vector
  \code{\link{createNewMap}} - creating de novo genetic map or chromosome assignment vector
  \code{\link{reduceChromosomesNumber}} - number of routines to reduce number of chromosomes of cross object
}

\keyword{manip}
