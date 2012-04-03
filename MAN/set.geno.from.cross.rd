\name{set.geno.from.cross}
\alias{set.geno.from.cross}


\title{Pull genotype from an object of class cross.}

\description{
  Pulling genotypes with a map from cross and putting into population object.
}

\usage{
set.geno.from.cross(cross,population,map=c("genetic","physical"))
	
}

\arguments{
\item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. If not supported, it will be created using data stored in population}
\item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details. }
  \item{map}{ 
  In which map ofd an population object shall map pulled from cross be stored:
  \itemize{
    \item{genetic}{ - genetic map - population$offspring$maps$genetic.}
    \item{physical}{ - physical map - population$offspring$maps$physical.}
  }
  }
}

\value{
  An object of class \code{\link{population}}. See \code{\link{create.population}} for details.
}

\details{
This function pull genotypes with a map from the cross object and puts them into provided population object. This is useful if the same cross is saturated and smoothed
(using \code{\link{smooth.geno}}) multiple times to fix the markers that were already smoothed.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
	data(yeastCross)
	yeastPopulation <- set.geno.from.cross(yeastCross,yeastPopulation)
}

\seealso{
  \itemize{
    \item{\code{\link{smooth.geno}}}{ - Remove genotyping errors by smoothing.}
    \item{\code{\link{reorganizeMarkersWithin}}}{ - Apply new ordering on the cross object usign ordering vector.}
    \item{\code{\link{assignChrToMarkers}}}{ - Create ordering vector from chromosome assignment vector.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or chromosome assignment vector.}
    \item{\code{\link{reduceChromosomesNumber}}}{ - Number of routines to reduce number of chromosomes of cross object.}
  }
}

\keyword{manip}
