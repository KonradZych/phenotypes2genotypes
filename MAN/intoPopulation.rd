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
	intoPopulation(population, dataObject, dataType=c("founders","offspring$phenotypes","offspring$genotypes","maps$genetic","maps$physical"),verbose=FALSE,debugMode=0)
}

\arguments{
 \item{population}{ object of class population, data should be put into}
 \item{dataObject}{ matrix of data to be put into ril object
 }
 \item{dataType}{ 
	specifies what kind of data dataObject contains:
	\itemize{
		\item{founders}{ - founders phenotype}
		\item{offspring$phenotypes }{ - offspring phenotype}
		\item{offspring$genotypes}{ - offspring genotype}
		\item{maps$genetic}{ - genetic map}
		\item{maps$physical}{ - physical map}
	}
 }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information }
}

\value{
  Object of class population. See \code{\link{createPopulation}} for description.
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
	population <- fakePopulation()
	offspring <- population$offspring$phenotypes
	founders <- population$founders$phenotypes
	population <- createPopulation(offspring,no.warn=TRUE)
	population <- intoPopulation(population,founders,"founders")
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{createPopulation}}
}

\keyword{manip}
