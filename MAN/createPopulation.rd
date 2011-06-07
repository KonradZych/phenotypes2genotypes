\name{createPopulation}
\alias{createPopulation}
\alias{population}


\title{Adding data to population object}

\description{
  Adding data to existing population object. Overwriting of existing data is permitted and will be executed with warning.
}

\usage{
	createPopulation(offspring_phenotypes=NULL, founders=NULL, offspring_genotypes=NULL, maps_genetic=NULL, maps_physical=NULL, no.warn=FALSE, verbose=FALSE, debugMode=0)
}

\arguments{
 \item{offspring_phenotypes}{ matrix containing offspring phenotype data (have to be supported, if not - function quits with error}
 \item{founders}{ matrix containing founders phenotype data (optional)}
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
	Konrad Zych \email{konrad.zych@uj.edu.pl}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
	Under tender patronage of: Danny Arends \email{Danny.Arends@gmail.com}
}

\examples{
	setwd(paste(.Library,"pheno2geno/data",sep="/"))
	founders <- read.table("founders_phenotypes.txt",sep="")
	offspring <- read.table("offspring_phenotypes.txt",sep="")
	population <- createPopulation(offspring)
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{intoPopulation}}
}

\keyword{manip}
