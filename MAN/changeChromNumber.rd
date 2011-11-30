\name{modify number of chromosomes}
\alias{reduceChromosomesNumber}
\alias{removeChromosomes}
\alias{removeTooSmallChromosomes}


\title{Change number of the chromosomes.}

\description{
  Methods to manually modify chromosomes number inside an object of class cross.
}

\usage{
	reduceChromosomesNumber(cross, numberOfChromosomes, verbose=FALSE)
	
	
	removeChromosomes(cross, chromosomesToBeRmv, verbose=FALSE)
	
	
	removeTooSmallChromosomes(cross, minNrOfMarkers, verbose=FALSE)
	
}

\arguments{
 \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
 \item{numberOfChromosomes}{ How many chromosomes should stay (remove all but 1:numberOfChromosomes).}
 \item{chromosomesToBeRmv}{ NAMES of chromosomes to be removed.}
 \item{minNrOfMarkers}{ Specify minimal number of markers chromosome is allowed to have (remove all that have less markers than that).}
  \item{verbose}{ Be verbose.}
}

\value{
  An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details.
}

\details{
There are three functions in pheno2geno enabling user to manually reduce number of resulting chromosomes. 

reduceChromosomesNumber
first of three functions is removing all but certain number of chromosomes. It depends only on ordering 
of chromosomes, not on naming. 


removeChromosomes
In opposite to previous one, this function operates only on names of chromosomes to be removed, not 
depending at all on the ordering.

removeTooSmallChromosomes
After using formLinkageGroups cross object often contains artifacts - linkage groups containing few 
markers. Those are probably really bad markers and should be removed. To do so, one could use 
removeTooSmallChromosomes function.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastCross)
	plot.rf(yeastCross, main="riself findBiomarkers example")
	cross_ <- reduceChromosomesNumber(yeastCross,5,verb=TRUE)
	plot.rf(cross_, main="Leaving only 5 chromosomes")
	cross_ <- removeChromosomes(yeastCross,1,verb=TRUE)
	plot.rf(cross_, main="Removing chromosome 1")
	cross_ <- removeTooSmallChromosomes(yeastCross,5,verb=TRUE)
	plot.rf(cross_, main="Leaving only chromosomes with more than 5 markers")
}

\seealso{
  \code{\link{reorganizeMarkersWithin}} - Apply new ordering on the cross object usign ordering vector.
  \code{\link{assignedChrToMarkers}} - Create ordering vector from chromosome assignment vector.
  \code{\link{cross.saturate}} - Saturate existing map.
  \code{\link{cross.denovo}} - Create de novo genetic map.
}

\keyword{manip}
