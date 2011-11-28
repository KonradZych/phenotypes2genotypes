\name{getBiomarkers}
\alias{getBiomarkers}

\title{Getting biomarkers.}

\description{
  Getting biomarkers from an object of class population.
}

\usage{
  getBiomarkers(population,pattern,verbose=FALSE)
}

\arguments{
 \item{population}{ an object of class \code{\link{population}}}
 \item{pattern}{ vector containg pattern to be matched in markers}
 \item{verbose}{ Be verbose}
}

\value{
  Matrix of all markers/ vector with markers matching best pattern.
}

\details{
	After running \code{\link{findBiomarkers}} function, biomarkers are stored inside \code{\link{population}} class object.
  To get them out, one can use getBiomarkers function. It will return matrix will all the markers or, if pattern is given,
  vector with marker matching best.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
  markers <- getBiomarkers(yeastPopulation,verbose=TRUE)
  bestMarker <- getBiomarkers(yeastPopulation,c(1,0,0,0,0,1,0,1,0,1,0,1,1,1,1,0,0),verbose=TRUE)
}

\seealso{
  \code{\link{readFiles}} - Loads genotype, phenotype, genetic map data files into R environment into a population object.
  \code{\link{createNewMap}} - Create de novo genetic map or vector showing how chromosomes should be assigned.
  \code{\link{saturateExistingMap}} - Saturate existing map.
  \code{\link{findDiffExpressed}} - Using Rank Product or student t-test analysis to select differentially expressed genes.
  \code{\link{findBiomarkers}} - Creating genotype markers  out of gene expression data.
}

\keyword{manip}
