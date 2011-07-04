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
  toGenotypes(population, genotype=c("simulated","real"), orderUsing=c("none","map_genetic","map_physical"), splitMethod=c("EM","mean"),treshold=0.01, overlapInd = 0, proportion = c(50,50), margin = 15, verbose=FALSE, debugMode=0,...)
}

\arguments{
 \item{population}{ Population type object, must contain parental phenotypic data.}
 \item{genotype}{ 
  Which genotypic matrix should be saved to file:
  \itemize{
    \item{simulated}{ - Genotype matrix from: \code{\link{toGenotypes}}}
    \item{real}{ - Original genotype matrix supplied by the user and read from file}
  }
  }
  \item{orderUsing}{ 
  which map should be used to order markers (Default - none, so markers are all put in 1 chromosome, with distance 1 cM between)
  \itemize{
    \item{map_genetic}{ - simulated data orderd using supplied genetic map}
    \item{map_physical}{ - simulated data orderd using supplied physical map}
  }
  }
 \item{splitMethod}{ Splitting markers using founders mean value or more sofisticated fitting of normal distributions by EM algoritm.}
 \item{treshold}{ If Rank Product pval for gene is lower that this value, we assume it is being diff. expressed.}
 \item{overlapInd}{ Number of individuals that are allowed in the overlap }
 \item{proportion}{ Proportion of individuals expected to carrying a certain genotype }
 \item{margin}{ Proportion is allowed to varry between this margin (2 sided) }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
  \item{...}{ Parameters passed to \code{\link[qtl]{formLinkageGroups}}. }
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
	cross <- toGenotypes(population,genotype="simulated",splitMethod="EM",proportion=c(25,50,25),orderUsing="map_genetic",treshold=0.1)
	plot.rf(cross, main="f2 toGenotypes example")

	population <- fakePopulation(type="bc")
	population <- findDiffExpressed(population)
	cross <- toGenotypes(population,genotype="simulated",splitMethod="EM",proportion=c(25,75),orderUsing="none",treshold=0.1)
	plot.rf(cross, main="bc toGenotypes example")

	population <- fakePopulation(type="riself")
	population <- findDiffExpressed(population)
	cross <- toGenotypes(population,genotype="simulated",splitMethod="mean",proportion=c(50,50),orderUsing="map_genetic",treshold=0.1)
	plot.rf(cross, main="riself toGenotypes example")
}

\seealso{
  \code{\link{readFiles}} - Loads genotype, phenotype, genetic map data files into R environment into a population object.
  \code{\link{findDiffExpressed}}
}

\keyword{manip}
