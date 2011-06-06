\name{intoPopulation}
\alias{intoPopulation}
\alias{intoPopulationSub.internal}
\alias{dataObject}
\alias{dataType}


\title{Adiing data to population object}

\description{
  Adding data to existing population object. Overwriting of existing data is permitted and will be executed with warning.
}

\usage{
	intoRil(population, dataObject, dataType=c("founders","offsprings"))
}

\arguments{
 \item{population}{ object of class population, data should be put into}
 \item{dataObject}{ matrix of data to be put into ril object}
 \item{dataType}{ what kind of data dataObject contains (parental phenotypic, children phenotypic)}
}

\value{
  Object of class population. See ?createPopulation for more details about structure.
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
	offsprings <- read.table("children_phenotypes.txt",sep="")
	population <- createPopulation(offsprings)
	population <- intoPopulation(population, founders, "founders")
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{createPopulation}}
}

\keyword{manip}
