\name{orderChromosomes}
\alias{orderChromosomes}


\title{Order chromosomes.}

\description{
  Order chromosomes of object of class cross using majority rule.
}

\usage{
	orderChromosomes(cross,method=c("majority","corelation"),map=c("genetic","physical"),verbose=FALSE)
}

\arguments{
 \item{cross}{ object of class cross}
 \item{method}{ 
  Which method should be used to compare (see details)
  \itemize{
    \item{majority}{ - majority}
    \item{corelation}{ - best corelation}
  }
  }
  \item{map}{ 
  Which map should be used for comparison:
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
	cross <- orderChromosomes(cross,method="majority",map="physical",verbose=TRUE)
	### not much happens, cause the are already ordered nicely
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{findDiffExpressed}}
  \code{\link{toGenotypes}}
}

\keyword{manip}
