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
 \item{ordering}{ ordering vector - for every marker in the cross object (names) specifies chromosome, this marker shall be moved to}
}

\value{
  an object of R/qtl class cross
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
  yeastCross <- reorganizeMarkersWithin(yeastCross, ordering)

}

\seealso{
  \code{\link{cross.denovo}} - creating de novo genetic map or chromosome assignment vector
  \code{\link{assignedChrToMarkers}} - create ordering vector from chromosome assignment vector
  \code{\link{cross.saturate}} - saturate existing map
}

\keyword{manip}

