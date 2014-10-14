\name{reorganizeMarkersWithin}
\alias{reorganizeMarkersWithin}

\title{Reorganize markers within cross object.}

\description{
  Reorder markers within cross object using supplied marker ordering vector.
}

\usage{
  reorganizeMarkersWithin(cross, ordering)
}

\arguments{
 \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
 \item{ordering}{ Ordering vector specifying for every marker in the cross object (by name) which new chromosome this marker will be moved to.}
}

\value{
  An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details.
}

\details{
  Functions reorders an object of class \code{cross} using the supplied marker ordering vector. This vector contains for each 
  marker the chromosome number that this marker will be moved to.
}

\author{
  Konrad Zych \email{k.zych@rug.nl}, Danny Arends \email{Danny.Arends@gmail.com}
  Maintainer: Konrad Zych \email{k.zych@rug.nl}
}

\examples{
  data(testCross)
  data(testPopulation)
  assignment <- cross.denovo(testPopulation,n.chr=5,verbose=TRUE,map="genetic",
  comparisonMethod=sumMajorityCorrelation, use.orderMarkers=FALSE,reOrder=FALSE)
  assignment #boring,but expected
  ordering <- assignChrToMarkers(assignment,testCross)
  testCross <- reorganizeMarkersWithin(testCross, ordering)

}

\seealso{
  \itemize{
    \item{\code{\link{cross.denovo}}}{ -  Creating a de novo genetic map or a chromosome assignment vector. }
    \item{\code{\link{cross.saturate}}}{ - Saturate an existing genetic map by using phenotype markers. }
    \item{\code{\link{assignChrToMarkers}}}{ - Create ordering vector from chromosome assignment vector. }
  }
}

\keyword{manip}

