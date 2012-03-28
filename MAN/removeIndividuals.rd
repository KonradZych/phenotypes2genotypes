\name{remove.individuals}
\alias{remove.individuals}


\title{Removing individuals.}

\description{
  Removing individuals from population object.
}

\usage{
	remove.individuals(population,individuals,verbose)
}

\arguments{
\item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details. }
 \item{individuals}{ individuals to be romved specified by their names}
  \item{verbose}{ Be verbose}
}

\value{
  An object of class \code{\link{population}}. See \code{\link{create.population}} for details.
}

\details{
  Function removes specified individuals from population object.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	### simulating data
	data(yeastPopulation)
	yeastPopulation <- remove.individuals(yeastPopulation,"V87",TRUE)

}

\seealso{
  \itemize{
    \item{\code{\link{read.population}}}{ - Load genotype, phenotype, genetic map data files into R environment into a population object.}
    \item{\code{\link{add.to.population}}}{ - Adding data to existing population object.}
  }
}

\keyword{manip}
