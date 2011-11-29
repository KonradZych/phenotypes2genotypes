\name{assignedChrToMarkers}
\alias{assignedChrToMarkers}


\title{Creating the ordering vector.}

\description{
 Create the ordering vector out of the chromosome assignment vector.
}

\usage{
	assignedChrToMarkers(assignment,cross)
}

\arguments{
 \item{assignment}{ chromosome assignment vector created using \link{cross.denovo} finction with reOrder = FALSE}
 \item{cross}{ an object of R/qtl class cross}
}

\value{
  Ordering vector, that can be used by \link{reorganizeMarkersWithin} function to reorder the cross object.
}

\details{
After using \link{cross.denovo} finction with reOrder = FALSE, chromosome assignment vector is created, showing how chromosomes
from created map shall be assigned to chromosomes from original map. This function uses this vector to create ordering vector, that can be used by 
\link{reorganizeMarkersWithin} function to reorder the cross object.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastCross)
	data(yeastPopulation)
	assignment <- cross.denovo(yeastPopulation,yeastCross,n.chr=16,verbose=TRUE,map="physical",comparisonMethod=sumMajorityCorrelation, use.orderMarkers=FALSE,reOrder=FALSE)
  assignment #boring,but expected
  ordering <- assignedChrToMarkers(assignment,yeastCross)
}

\seealso{
  \code{\link{reorganizeMarkersWithin}} - apply new ordering on the cross object usign ordering vector
  \code{\link{cross.saturate}} - saturate existing map
  \code{\link{cross.denovo}} - create de novo genetic map or vector showing how chromosomes should be assigned
}

\keyword{manip}
