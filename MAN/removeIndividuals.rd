\name{removeIndividuals}
\alias{removeIndividuals}


\title{Removing individuals from population object.}

\description{
  Removing specified individuals from population object.
}

\usage{
	removeIndividuals(population,individuals)
}

\arguments{
 \item{population}{ object of class population}
 \item{individuals}{ individuals to be romved specified by their names}
}

\value{
  Cross type object.
}

\details{
  TODO
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
	Under tender patronage of: Danny Arends \email{Danny.Arends@gmail.com}
}

\examples{
	### simulating data
	population <- fakePopulation()
	population <- removeIndividuals(population,"RIL_1")

}

\seealso{
  \code{\link{readFiles}}
  \code{\link{preprocessData}}
  \code{\link{toGenotypes}}
}

\keyword{manip}
