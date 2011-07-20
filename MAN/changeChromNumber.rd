\name{modify number of chromosomes}
\alias{reduceChromosomesNumber}
\alias{removeChromosomes}
\alias{removeTooSmallChromosomes}


\title{Methods to modify chromosomes number.}

\description{
  Methods to modify chromosomes number of cross class object.
}

\usage{
	reduceChromosomesNumber(cross, numberOfChromosomes, verbose=FALSE)
	
	
	removeChromosomes(cross, chromosomesToBeRmv, verbose=FALSE)
	
	
	removeTooSmallChromosomes(cross, minNrOfMarkers, verbose=FALSE)
	
}

\arguments{
 \item{cross}{ object of class cross}
 \item{numberOfChromosomes}{ how many chromosomes should stay (remove all but 1:numberOfChromosomes)}
 \item{chromosomesToBeRmv}{ explicitly provide functions with NAMES of chromosomes to be removed}
 \item{minNrOfMarkers}{ specify minimal number of markers chromosome is allowed to have (remove all that have less markers than that)}
  \item{verbose}{ be verbose}
}

\value{
  Cross type object.
}

\details{
  TODO
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	set.seed(102)
	population <- fakePopulation(type="riself",n.markers=2000)
	population <- findDiffExpressed(population)
	cross <- toGenotypes(population,genotype="simulated",proportion=c(50,50),orderUsing="map_genetic",treshold=0.01)
	plot.rf(cross, main="riself toGenotypes example")
	cross_ <- reduceChromosomesNumber(cross,5,verb=TRUE)
	plot.rf(cross_, main="Leaving only 5 chromosomes")
	cross_ <- removeChromosomes(cross,1,verb=TRUE)
	plot.rf(cross_, main="Removing chromosome 1")
	cross_ <- removeTooSmallChromosomes(cross,5,verb=TRUE)
	plot.rf(cross_, main="Leaving only chromosomes with more than 5 markers")


}

\seealso{
  \code{\link{readFiles}}
  \code{\link{findDiffExpressed}}
  \code{\link{toGenotypes}}
}

\keyword{manip}
