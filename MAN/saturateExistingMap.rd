\name{Saturate existing map}
\alias{cross.saturate}


\title{Saturate existing map.}

\description{
  Saturating existing map.
}

\usage{
	cross.saturate(population, cross, map=c("genetic","physical"), corSDTreshold=3, use.orderMarkers=FALSE, verbose=FALSE, debugMode=0)
	
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
 \item{corSDTreshold}{ specifies threshold for selecting/rejecting markers (see \link{markerPlacementPlot} )}
 \item{use.orderMarkers}{should markers on the newly created map be ordered using R/qtl orderMarkers funtion}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information }
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
	cross <- cross.saturate(yeastPopulation,map="physical",verbose=TRUE,debugMode=2)
}

\seealso{
  \code{\link{reorganizeMarkersWithin}} - apply new ordering on the cross object usign ordering vector
  \code{\link{assignedChrToMarkers}} - create ordering vector from chromosome assignment vector
  \code{\link{cross.denovo}} - creating de novo genetic map or chromosome assignment vector
  \code{\link{reduceChromosomesNumber}} - number of routines to reduce number of chromosomes of cross object
  \code{\link{markerPlacementPlot}} - plot showing how many markers will be selected for map saturation with different thresholds
}

\keyword{manip}
