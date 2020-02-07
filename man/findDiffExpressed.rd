\name{find.diff.expressed}
\alias{find.diff.expressed}
\alias{groupLabels}

\title{Finding differentially expressed genes.}

\description{
  Using Rank Product or student t-test analysis to select differentially expressed genes.
}

\usage{
  find.diff.expressed(population, verbose=FALSE, debugMode=0,...)
}

\arguments{
 \item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details. }
 \item{verbose}{ Be verbose.}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information.}
 \item{...}{ Additional arguments passed to RP function.}
}

\value{
  Object of class \code{\link{population}}.
}

\details{
  This function finds probes differentially expressed between founders using either student t.test
}

\author{
	Konrad Zych \email{k.zych@rug.nl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{k.zych@rug.nl}
}

\examples{
	data(testPopulation)
	testPopulation <- find.diff.expressed(testPopulation)
}

\seealso{
  \itemize{
    \item{\code{\link{read.population}}}{ - Load genotype, phenotype, genetic map data files into R environment into a population object.}
    \item{\code{\link{generate.biomarkers}}}{ - Creating genotypes from children phenotypes.}
  }
}

\keyword{manip}
