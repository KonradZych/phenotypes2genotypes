\name{cleanMap}
\alias{cleanMap}


\title{Removing wrong markers.}

\description{
  Removing markers that cause reco map to be too long.
}

\usage{
	cleanMap(cross, difPercentage, minChrLenght, verbose=FALSE, debugMode=0)
}

\arguments{
 \item{cross}{ R/qtl cross type object.}
 \item{difPercentage}{ If removing marker will make map shorter by this percentage of its length, then it will be dropped.}
 \item{minChrLenght}{ chromosomes shorter than that won't be processed}
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
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
	### simulating data
	population <- fakePopulation()
	cross <- toGenotypes(population,genotype="real",orderUsing="map_genetic")
	### dangerous
	cross <- cleanMap(cross,10,70,verbose=TRUE, debugMode=1)
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{preprocessData}}
  \code{\link{toGenotypes}}
}

\keyword{manip}
