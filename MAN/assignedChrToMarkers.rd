\name{assignChrToMarkers}
\alias{assignChrToMarkers}

\title{Create an ordering vector to assign markers to chromosomes}

\description{
 Creates an ordering vector out of the chromosome assignment vector.
}

\usage{
	assignChrToMarkers(assignment, cross)
}

\arguments{
 \item{assignment}{ Chromosome assignment vector created using \link{cross.denovo} with reOrder = FALSE}
 \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
}

\value{
  Ordering vector, that can be used by \link{reorganizeMarkersWithin} function to reorder the cross object.
}

\details{
When using \link{cross.denovo} function with reOrder = FALSE, an chromosome assignment vector is created. 
This vector shows how chromosomes from the created map are assigned to chromosomes from the original map. 
This function transforms this vector and creates an ordering vector, that can be used by the
\code{\link{reorganizeMarkersWithin}} function to reorder markers inside the cross object.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastCross)
	data(yeastPopulation)
	assignment <- cross.denovo(yeastPopulation,n.chr=16,verbose=TRUE,map="physical",comparisonMethod=sumMajorityCorrelation, use.orderMarkers=FALSE,reOrder=FALSE, cross=yeastCross)
  assignment #boring,but expected
  ordering <- assignChrToMarkers(assignment,yeastCross)
}

\seealso{
  \code{\link{reorganizeMarkersWithin}} - Apply new ordering on the cross object usign ordering vector.
  \code{\link{cross.saturate}} - Saturate existing map.
  \code{\link{cross.denovo}} - Create de novo genetic map or vector showing how chromosomes should be assigned.
}

\keyword{manip}
