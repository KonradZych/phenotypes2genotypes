\name{rearrangeMarkers}
\alias{rearrangeMarkers}


\title{Rearrangeing markers.}

\description{
  Rearrangeing markers inside an object of class cross using correlation.
}

\usage{
	rearrangeMarkers(cross,map=c("genetic","physical"),corTreshold=0.6,addMarkers=FALSE,verbose=FALSE)
}

\arguments{
 \item{cross}{ object of class cross}
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
  an object of R/qtl class cross
}

\details{
  TODO
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	### simulating data
	population <- fakePopulation(type="riself",n.markers=1000)
	population <- findDiffExpressed(population)
	cross <- toGenotypes(population,genotype="simulated",proportion=c(50,50),orderUsing="map_genetic",treshold=0.01)
	plot.rf(cross)
	#cross <- rearrangeMarkers(cross,map="physical",verbose=TRUE)
	#plot.rf(cross)
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{findDiffExpressed}}
  \code{\link{toGenotypes}}
}

\keyword{manip}