\name{reorganizeMarkersWithin}
\alias{reorganizeMarkersWithin}


\title{Reorganize markers within cross object.}

\description{
  Reorders markers withion cross object using supplied ordering vector.
}

\usage{
	reorganizeMarkersWithin(cross, ordering)
}

\arguments{
 \item{cross}{ object of class population}
 \item{ordering}{ number of chromosomes expected}
}

\value{
  an object of R/qtl class cross
}


\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
  data(yeastCross)
	assignment <- createNewMap(yeastPopulation,16,verbose=TRUE,map="physical",reOrder=FALSE,comparisonMethod=sumMajorityCorrelation, use.orderMarkers=FALSE)
  ordering <- assignedChrToMarkers(assignment,yeastCross)
  yeastCross <- reorganizeMarkersWithin(yeastCross, ordering)

}

\seealso{
  \code{\link{createNewMap}} - creating de novo genetic map or chromosome assignment vector
  \code{\link{assignedChrToMarkers}} - create ordering vector from chromosome assignment vector
  \code{\link{saturateExistingMap}} - saturate existing map
}

\keyword{manip}

