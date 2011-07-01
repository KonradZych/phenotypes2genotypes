\name{orderChromosomes}
\alias{orderChromosomes}


\title{Order chromosomes.}

\description{
  Order chromosomes of object of class cross using majority rule.
}

\usage{
	orderChromosomes(cross,map=c("genetic","physical"),verbose=FALSE)
}

\arguments{
 \item{cross}{ object of class cross}
 \item{map}{ 
  Which genotypic matrix should be saved to file:
  \itemize{
    \item{genetic}{ - genetic map from cross$maps$genetic}
    \item{physical}{ - physical map from cross$maps$physical}
  }
  }
 \item{verbose}{ Be verbose}
}

\value{
  object of class cross
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
	cross <- orderChromosomes(cross,"physical",verbose=TRUE)
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{findDiffExpressed}}
  \code{\link{toGenotypes}}
}

\keyword{manip}
