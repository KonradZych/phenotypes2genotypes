\name{Create new map}
\alias{createNewMap}


\title{Creating de novo genetic map.}

\description{
  Creating de novo genetic map.
}

\usage{
	enrichExistingMap(population,cross,map=c("genetic","physical"),corTreshold=0.6,verbose=FALSE,debugMode=0)
	
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
This function enriches existing map (stored in population object) with markers derived form gene expression data
(provided inside cross or population. Correlation matrix between those two is made and based on it, new markers sre
being placed on the map.
}
\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	#population <- fakePopulation()
	#population <- findDiffExpressed(population)
	#population <- toGenotypes(population)
	#cross <-enrichExistingMap(population)
}

\seealso{
  \code{\link{orderChromosomes}} - ordering chromosomes of an object of class cross using majority rule
  \code{\link{rearrangeMarkers}} - rearrangeing markers inside an object of class cross using correlation
  \code{\link{reduceChromosomesNumber}} - removing all but certain number of chromosomes from an object of class cross
}

\keyword{manip}
