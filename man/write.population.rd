\name{write.population}
\alias{write.population}

\title{Writes a population object to file.}

\description{
  Writes a population object to file, for easy loading of intermediate data later.
}

\usage{
	write.population(population, offspring = "offspring", founders = "founders", map = "map", verbose = FALSE, debugMode = 0)
}

\arguments{
 \item{population}{ An object of class \code{population}. See \code{\link{create.population}} for details. }
 \item{offspring}{ Core used to specify names of children phenotypic ("core_phenotypes.txt") genotypic ("core_genotypes.txt") and annotations ("core_annotations.txt") files.}
 \item{founders}{ Core used to specify names of parental phenotypic ("core_phenotypes.txt") file. }
 \item{map}{ Core used to specify names of genetic ("map_genetic.txt") and physical ("map_physical.txt") map files. }
 \item{verbose}{ Be verbose. }
 \item{debugMode}{ 1: Print out checks, 2: print additional time information }
}

\value{
  None.
}

\details{
  This function writes an object of class \code{population} into a file.
}

\author{
	Konrad Zych \email{k.zych@rug.nl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{k.zych@rug.nl}
}

\examples{
	\dontrun{
		population <- fake.population()
		write.population(population,verbose=TRUE)
	}
}

\seealso{
  \itemize{
    \item{\code{\link{add.to.population}}}{ - Adding data to existing population object.}
    \item{\code{\link{create.population}}}{ - Create new object of class population.}
    \item{\code{\link{read.population}}}{ - Create new object of class population.}
  }
}

\keyword{manip}
