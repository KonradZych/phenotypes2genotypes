\name{modify number of chromosomes}
\alias{reduceChromosomesNumber}
\alias{removeChromosomes}
\alias{removeTooSmallChromosomes}

\title{Change the number of chromosomes in a cross object}

\description{
  Methods to manually modify the number of chromosomes inside an cross object
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
There are three functions in pheno2geno to allow the user to manually reduce number of resulting chromosomes. 

\emph{reduceChromosomesNumber}
This functions removes all chromosomes from the cross object excluding chromosome 1 to numberOfChromosomes. It depends on the ordering 
of chromosomes inside the cross object (which is based on the length of the chromosomes). 

\emph{ removeChromosomes}
This function removes chromosomes from the cross object by name. Because of this it does not depend on the 
ordering of the chromosomes inside the cross object.

\emph{ removeTooSmallChromosomes}
This function is used to clean a cross object after using \code{\link{formLinkageGroups}}. FormLinkageGroups can 
introduce small chromosomes as artifacts. These linkage groups consist of only a few markers with poor quality 
and should be removed from the cross object.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastCross)
	plot.rf(yeastCross, main="riself generate.biomarkers example")
	cross_ <- reduceChromosomesNumber(yeastCross,5,verb=TRUE)
	plot.rf(cross_, main="Leaving only 5 chromosomes")
	cross_ <- removeChromosomes(yeastCross,1,verb=TRUE)
	plot.rf(cross_, main="Removing chromosome 1")
	cross_ <- removeTooSmallChromosomes(yeastCross,5,verb=TRUE)
	plot.rf(cross_, main="Leaving only chromosomes with more than 5 markers")
}

\seealso{
  \code{\link{reorganizeMarkersWithin}} - Apply new ordering on the cross object usign ordering vector.
  \code{\link{assignChrToMarkers}} - Create ordering vector from chromosome assignment vector.
  \code{\link{cross.saturate}} - Saturate existing map.
  \code{\link{cross.denovo}} - Create de novo genetic map.
}

\keyword{manip}
