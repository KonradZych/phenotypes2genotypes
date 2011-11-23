\name{toGenotypes}
\alias{toGenotypes}
\alias{use}
\alias{treshold}
\alias{overlapInd}
\alias{proportion}
\alias{margin}
\alias{splitMethod}
\alias{numberOfChromosomes}

\title{Creating genotypes from children phenotypes}

\description{
  Creating genotype markers  out of gene expression data.
}

\usage{
  toGenotypes(population, treshold=0.05, overlapInd = 0, proportion = c(50,50), margin = 15, verbose=FALSE, debugMode=0)
}

\arguments{
 \item{population}{ Population type object, must contain parental phenotypic data.}
 \item{treshold}{ If Rank Product pval for gene is lower that this value, we assume it is being diff. expressed.}
 \item{overlapInd}{ Number of individuals that are allowed in the overlap }
 \item{proportion}{ Proportion of individuals expected to carrying a certain genotype }
 \item{margin}{ Proportion is allowed to varry between this margin (2 sided) }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information }
}

\value{
  Cross type object. 
}

\details{
	This function, using Mixture Models splits offspring gene expression data into genotype markers, based on founders gene expression data.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	#Example for F2 population
	set.seed(102)
	population <- fakePopulation(type="f2")
	population <- findDiffExpressed(population)
	population <- toGenotypes(population,proportion=c(25,50,25),treshold=0.01)
	\dontrun{
	#Example for BC population
	set.seed(102)
	population <- fakePopulation(type="bc")
	population <- findDiffExpressed(population)
	population <- toGenotypes(population,proportion=c(25,75),treshold=0.01)

	#Example for BC population
	set.seed(102)
	population <- fakePopulation(type="riself")
	population <- findDiffExpressed(population)
	population <- toGenotypes(population,proportion=c(50,50),treshold=0.01)
	}
}

\seealso{
  \code{\link{readFiles}} - loads genotype, phenotype, genetic map data files into R environment into a population object.
  \code{\link{createNewMap}} - create de novo genetic map or vector showing how chromosomes should be assigned
  \code{\link{saturateExistingMap}} - saturate existing map
  \code{\link{findDiffExpressed}} - using Rank Product or student t-test analysis to select differentially expressed genes
}

\keyword{manip}
