\name{reorganizeMarkersWithin}
\alias{reorganizeMarkersWithin}


\title{Reorganize markers within cross object.}

\description{
  Reorders markers withion cross object using supplied ordering vector.
}

\usage{
	reorganizeMarkersWithin(cross, ordering)
}

\arguments{
 \item{cross}{ object of class population}
 \item{ordering}{ number of chromosomes expected}
}

\value{
  an object of R/qtl class cross
}

\details{
  }

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	### simulating data
	population <- fakePopulation(type="riself")
	population <- findDiffExpressed(population)
	population <- toGenotypes(population,genotype="simulated",proportion=c(50,50),orderUsing="map_genetic",treshold=0.01)
	cross <- assignChromosomes(population,10,map="physical",verbose=TRUE)

}

\seealso{
  \code{\link{readFiles}}
  \code{\link{findDiffExpressed}}
  \code{\link{toGenotypes}}
}

\keyword{manip}

