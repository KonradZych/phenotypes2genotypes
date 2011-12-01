\name{read.population}
\alias{read.population}

\title{Readinf population object from file.}

\description{
  Loading object of class population from text file.
}

\usage{
	read.population(filename="population.txt",verbose=FALSE)
}

\arguments{
 \item{filename}{ Name of the population file.}
 \item{verbose}{ Be verbose.}
}

\value{
  An object of class \code{\link{population}}. See \code{\link{createPopulation}} for details.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	\dontrun{
	read.population("population.txt")
	}
}

\seealso{
	\itemize{
    \item{\code{\link{readFiles}}}{ - Load genotype, phenotype, genetic map data files into R environment into a population object.}
    \item{\code{\link{createPopulation}}}{ - Create new object of class population.}
    \item{\code{\link{intoPopulation}}}{ - Adding data to existing population object.}
  }
}

\keyword{manip}
