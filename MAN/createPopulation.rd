\name{createPopulation}
\alias{createPopulation}
\alias{intoPopulationSub.internal}
\alias{dataObject}
\alias{dataType}


\title{Adiing data to population object}

\description{
  Adding data to existing population object. Overwriting of existing data is permitted and will be executed with warning.
}

\usage{
	createPopulation(offspring_phenotypes=NULL, founders=NULL, offspring_genotypes=NULL, maps_genetic=NULL, maps_physical=NULL, verbose=FALSE, debugMode=0)
}

\arguments{
 \item{offspring_phenotypes}{ matrix containing offspring phenotype data (have to be supported, if not - function quits with error}
 \item{founders}{ matrix containing founders phenotype data (optional)}
 \item{offspring_genotypes}{ matrix containing offspring genotype data (optional)}
 \item{maps_genetic}{ matrix containing genetic map (optional)}
 \item{maps_physical}{ matrix containing physical map (optional)}
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
	Konrad Zych \email{konrad.zych@uj.edu.pl}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
	Under tender patronage of: Danny Arends \email{Danny.Arends@gmail.com}
}

\examples{
	setwd(paste(.Library,"pheno2geno/data",sep="/"))
	founders <- read.table("parental_phenotypes.txt",sep="")
	offspring <- read.table("children_phenotypes.txt",sep="")
	population <- createPopulation(offspring)
	population <- intoPopulation(population, founders, "founders")
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{intoPopulation}}
}

\keyword{manip}
