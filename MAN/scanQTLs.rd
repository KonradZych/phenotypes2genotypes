\name{scanQTLs}
\alias{scanQTLs}

\title{Scan qtls}

\description{
  Scanning population data for qtls for use in cross.saturate function.
}

\usage{
  scanQTLs(population,map=c("genetic","physical"),step=0.1,verbose=FALSE)
}

\arguments{
  \item{population}{ An object of class \code{\link{population}}. See \code{\link{createPopulation}} for details. }
   \item{map}{ Which map (from ones stored in population$maps) should be used fo assigning chromosomes on the created map}
 \item{step}{ Maximum distance (in cM) between positions at which the genotype probabilities are calculated, though for step = 0, probabilities are calculated only at the marker locations. See \code{\link{calc.genoprob} for more information} .}
 \item{verbose}{ Be verbose.}
}

\value{
  An object of class \code{\link{population}}. See \code{\link{createPopulation}} for details.
}

\details{
	This function takes care about qtl scan that is used by \code{\link{cross.saturate}} function. It was made separated function, since process itself takes a long time
  and before running \code{\link{cross.saturate}} function one should run \code{\link{markerPlacementPlot}} to assess the optimal threshold.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
	yeastPopulation <- scanQTLs(yeastPopulation,verbose=TRUE,map="physical",step=2)
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
