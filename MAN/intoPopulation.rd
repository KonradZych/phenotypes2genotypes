\name{intoPopulation}
\alias{intoPopulation}
\alias{intoPopulationSub.internal}
\alias{dataObject}
\alias{dataType}


\title{Adding data to population object}

\description{
  Adding data to existing population object. Overwriting of existing data is permitted and will be executed with warning.
}

\usage{
	intoPopulation(population, dataObject, dataType=c("founders","offspring$phenotypes","founders$group","offspring$genotypes","maps$genetic","maps$physical"),verbose=FALSE,debugMode=0)
}

\arguments{
 \item{population}{ An object of class \code{\link{population}}. See \code{\link{createPopulation}} for details. }
 \item{dataObject}{ A matrix of data to be put into ril object or a list of matrices (then dataType should be a list as well).}
 \item{dataType}{ 
	Specifies what kind of data dataObject contains:
	\itemize{
		\item{founders}{ - Founders phenotype.}
		\item{offspring$phenotypes }{ - Offspring phenotype.}
		\item{founders$group }{ - Specifying groups in founders phenotypes.}
		\item{offspring$genotypes}{ - Offspring genotype.}
		\item{maps$genetic}{ - Genetic map.}
		\item{maps$physical}{ - Physical map.}
	}
 }
 \item{verbose}{ Be verbose.}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information.}
}

\value{
  An object of class \code{\link{population}}. See \code{\link{createPopulation}} for details.
}

\details{
  This function inputs data into existing population object. It can input single matrix or list of matrices.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	population <- fakePopulation()
	offspring <- population$offspring$phenotypes
	founders <- population$founders$phenotypes
	founders_groups <- population$founders$groups
	maps_genetic <- population$maps$genetic
	population <- createPopulation(offspring,founders,founders_groups)
	population <- intoPopulation(population,maps_genetic,"maps$genetic")
}

\seealso{
  \itemize{
    \item{\code{\link{readFiles}}}{ - Loads genotype, phenotype, genetic map data files into R environment into a population object.}
    \item \code{\link{createPopulation}}}{ - Create object of class population from data already in R environment.}
    \item \code{\link{fakePopulation}}}{ - Simulate basic population object for use in examples.}
  }
}

\keyword{manip}
