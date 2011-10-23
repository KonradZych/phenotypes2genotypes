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
  Creating genotypes from children phenotypes using parental data and saving cross object.
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
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
}

\value{
  Cross type object. 
}

\details{
	TODO
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	set.seed(102)
	population <- fakePopulation(type="f2")
	population <- findDiffExpressed(population)
	cross <- toGenotypes(population,genotype="simulated",proportion=c(25,50,25),orderUsing="map_genetic",treshold=0.01)
	plot.rf(cross, main="f2 toGenotypes example")

	population <- fakePopulation(type="bc")
	population <- findDiffExpressed(population)
	cross <- toGenotypes(population,genotype="simulated",proportion=c(25,75),orderUsing="none",treshold=0.01)
	plot.rf(cross, main="bc toGenotypes example")

	population <- fakePopulation(type="riself")
	population <- findDiffExpressed(population)
	cross <- toGenotypes(population,genotype="simulated",proportion=c(50,50),orderUsing="map_genetic",treshold=0.01)
	plot.rf(cross, main="riself toGenotypes example")
}

\seealso{
  \code{\link{readFiles}} - Loads genotype, phenotype, genetic map data files into R environment into a population object.
  \code{\link{findDiffExpressed}}
}

\keyword{manip}
