\name{smoothGeno}
\alias{smoothGeno}


\title{Smooth genotype.}

\description{
 Remove genotyping errors by smoothing.
}

\usage{
	smoothGeno(cross,windowSize=1,chr,verbose=FALSE)
}

\arguments{
 \item{cross}{ an object of R/qtl class cross}
 \item{windowSize}{ specifies number of markers that are processed at once}
 \item{chr}{ specifies subset of chromosomes to be processed}
  \item{verbose}{ be verbose}
}

\value{
 An object of R/qtl class cross.
}

\details{
This function is usefull while creating new map (also a bit while saturating existing map). It is taking a number of markers at the time (specified by user) and checking values 
of markers bordering those selected. If border values at both sides are the same and they are different from value of selected markers, then selected markers are given value of
border ones. Then the windows slides by one marker and another group is chacked.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastCross)
  yeastCross_smooth <- smoothGeno(yeastCross,3)
  geno.image(yeastCross)
	geno.image(yeastCross_smooth)
  
}

\seealso{
  \code{\link{saturateExistingMap}} - saturate existing map
  \code{\link{createNewMap}} - create de novo genetic map or vector showing how chromosomes should be assigned
}

\keyword{manip}
