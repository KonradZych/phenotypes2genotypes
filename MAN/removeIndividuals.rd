\name{removeIndividuals}
\alias{removeIndividuals}


\title{Removing individuals from population object.}

\description{
  Removing specified individuals from population object.
}

\usage{
	removeIndividuals(population,individuals,verbose)
}

\arguments{
 \item{population}{ object of class population}
 \item{individuals}{ individuals to be romved specified by their names}
  \item{verbose}{ Be verbose}
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
	population <- removeIndividuals(population,"RIL_1",TRUE)

}

\seealso{
  \code{\link{readFiles}}
  \code{\link{findDiffExpressed}}
  \code{\link{toGenotypes}}
}

\keyword{manip}
