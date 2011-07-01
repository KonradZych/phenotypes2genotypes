\name{createPopulation}
\alias{createPopulation}
\alias{population}


\title{Create new object of class population.}

\description{
  Create new object of class population. If object exists, will be overwritten.
}

\usage{
	createPopulation(offspring_phenotypes, founders, founders_groups, offspring_genotypes, maps_genetic, maps_physical, no.warn=FALSE, verbose=FALSE,debugMode=0)
}

\arguments{
 \item{offspring_phenotypes}{ matrix containing offspring phenotype data (have to be supported, if not - function quits with error)}
 \item{founders}{ matrix containing founders phenotype data (have to be supported, if not - function quits with error)}
  \item{founders_groups}{ specify groups im founders data (have to be supported, if not - function quits with error)}
 \item{offspring_genotypes}{ matrix containing offspring genotype data (optional)}
 \item{maps_genetic}{ matrix containing genetic map (optional)}
 \item{maps_physical}{ matrix containing physical map (optional)}
 \item{no.warn}{ if TRUE, no warnings will be produced}
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information }
}

\value{
  Object of class population. TODO: description
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
	offspring <- population$offspring$phenotypes
	founders <- population$founders$phenotypes
	founders_groups <- population$founders$groups
	population <- createPopulation(offspring,founders,founders_groups)
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{intoPopulation}}
}

\keyword{manip}
