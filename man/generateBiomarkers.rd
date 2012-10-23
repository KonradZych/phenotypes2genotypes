\name{generate.biomarkers}
\alias{generate.biomarkers}
\alias{use}
\alias{threshold}
\alias{overlapInd}
\alias{proportion}
\alias{margin}
\alias{splitMethod}
\alias{numberOfChromosomes}

\title{Generate discrete biomarkers from the continuous phenotypes}

\description{
  Creating genotype markers  out of gene expression data.
}

\usage{
  generate.biomarkers(population, threshold=0.05, overlapInd = 10, proportion = c(50,50), margin = 15, verbose=FALSE, debugMode=0)
}

\arguments{
  \item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details. }
  \item{threshold}{ If the pvalue for differential expression of this phenotype (see \code{\link{find.diff.expressed}}) is 
  lower that the set threshold, the phenotype is kept in the analysis as being differentially expressed.}
  \item{overlapInd}{ The number of individuals that are allowed in the overlap (undecided region) when assigning genotype encodings.}
  \item{proportion}{ The expected proportion of individuals expected to carrying a certain genotype (e.g. c(50,50) in a recombinant inbred line).}
  \item{margin}{ This specifies how much deviation from the expected proportion is allowed (2 sided). }
  \item{verbose}{ Be verbose. }
  \item{debugMode}{ Either use 1 or 2, this will modify the amount of information returned to the user. 1) Print out checks, 2) Print additional time information.}
}

\value{
  An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details 
}

\details{
	This function, using the results from mixture modeling splits the continuous offspring phenotype data into discrete genotype markers, infering the 
  direction from the founders expression data.
}

\author{
	Konrad Zych \email{k.zych@rug.nl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{k.zych@rug.nl}
}

\examples{
	#Example for F2 population
	set.seed(102)
	population <- fake.population(type="f2")
	population <- find.diff.expressed(population)
	population <- generate.biomarkers(population,proportion=c(25,50,25),threshold=0.01)
	\dontrun{
	#Example for BC population
	set.seed(102)
	population <- fake.population(type="bc")
	population <- find.diff.expressed(population)
	population <- generate.biomarkers(population,proportion=c(25,75),threshold=0.01)

	#Example for BC population
	set.seed(102)
	population <- fake.population(type="riself")
	population <- find.diff.expressed(population)
	population <- generate.biomarkers(population,proportion=c(50,50),threshold=0.01)
	}
}

\seealso{
  \itemize{
    \item{\code{\link{read.population}}}{ - Load genotype, phenotype, genetic map data files into R environment into a population object.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or vector showing how chromosomes should be assigned.}
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{find.diff.expressed}}}{ - Using Rank Product or student t-test analysis to select differentially expressed genes.}
  }
}

\keyword{manip}
