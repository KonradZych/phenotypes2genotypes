\name{findBiomarkers}
\alias{findBiomarkers}
\alias{use}
\alias{threshold}
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
  findBiomarkers(population, threshold=0.05, overlapInd = 0, proportion = c(50,50), margin = 15, verbose=FALSE, debugMode=0)
}

\arguments{
  \item{population}{ An object of class \code{\link{population}}. See \code{\link{createPopulation}} for details. }
 \item{threshold}{ If pval for gene (see \code{\link{findDiffExpressed}}) is lower that this value, we assume it is being diff. expressed.}
 \item{overlapInd}{ Number of individuals that are allowed to overlap between genotypes.}
 \item{proportion}{ Proportion of individuals expected to carrying a certain genotype.}
 \item{margin}{ Proportion is allowed to varry between this margin (2 sided).}
 \item{verbose}{ Be verbose.}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information.}
}

\value{
  An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details 
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
	population <- findBiomarkers(population,proportion=c(25,50,25),treshold=0.01)
	\dontrun{
	#Example for BC population
	set.seed(102)
	population <- fakePopulation(type="bc")
	population <- findDiffExpressed(population)
	population <- findBiomarkers(population,proportion=c(25,75),treshold=0.01)

	#Example for BC population
	set.seed(102)
	population <- fakePopulation(type="riself")
	population <- findDiffExpressed(population)
	population <- findBiomarkers(population,proportion=c(50,50),treshold=0.01)
	}
}

\seealso{
  \itemize{
    \item{\code{\link{readFiles}}}{ - Load genotype, phenotype, genetic map data files into R environment into a population object.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or vector showing how chromosomes should be assigned.}
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{findDiffExpressed}}}{ - Using Rank Product or student t-test analysis to select differentially expressed genes.}
  }
}

\keyword{manip}
