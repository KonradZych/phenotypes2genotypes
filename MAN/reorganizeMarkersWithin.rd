\name{reorganizeMarkersWithin}
\alias{reorganizeMarkersWithin}


\title{Reorganize markers within cross object.}

\description{
  Reorder markers within cross object using supplied ordering vector.
}

\usage{
	reorganizeMarkersWithin(cross, ordering)
}

\arguments{
 \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
 \item{ordering}{ Ordering vector - for every marker in the cross object (names) specifies chromosome, this marker shall be moved to.}
}

\value{
  An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details.
}

\details{
  Functions reorders an object of class \code{cross} using supplied vector containing, for each of the markers, chromosome number this marker shall be moved to.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastCross)
	data(yeastPopulation)
	assignment <- cross.denovo(yeastPopulation,n.chr=16,verbose=TRUE,map="physical",comparisonMethod=sumMajorityCorrelation, cross = yeastCross, use.orderMarkers=FALSE,reOrder=FALSE)
  assignment #boring,but expected
  ordering <- assignChrToMarkers(assignment,yeastCross)
  yeastCross <- reorganizeMarkersWithin(yeastCross, ordering)

}

\seealso{
  \itemize{
    \item{\code{\link{cross.denovo}}}{ -  Creating de novo genetic map or chromosome assignment vector.}
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{assignChrToMarkers}}}{ - Create ordering vector from chromosome assignment vector.}
  }
}

\keyword{manip}

