\name{remove.individuals}
\alias{remove.individuals}


\title{Remove individuals from a population object.}

\description{
  Remove individuals from a population object.
}

\usage{
	remove.individuals(population, individuals,verbose)
}

\arguments{
  \item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details. }
  \item{individuals}{ individuals to be removed specified by name }
  \item{verbose}{ Be verbose. }
}

\value{
  An object of class \code{\link{population}}. See \code{\link{create.population}} for details.
}

\details{
  This function removes the specified individuals from a population object.
}

\author{
	Konrad Zych \email{k.zych@rug.nl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{k.zych@rug.nl}
}

\examples{
	### simulating data
	data(yeastPopulation)
	yeastPopulation <- remove.individuals(yeastPopulation, "V87", TRUE)
}

\seealso{
  \itemize{
    \item{\code{\link{read.population}}}{ - Load genotype, phenotype, genetic map data files into R environment into a population object.}
    \item{\code{\link{add.to.population}}}{ - Adding data to existing population object.}
  }
}

\keyword{manip}
