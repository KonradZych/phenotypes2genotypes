\name{cross.saturate}
\alias{cross.saturate}


\title{Saturate existing map.}

\description{
  Saturating existing map usign markers derived from geene expression data.
}

\usage{
cross.saturate(population, cross, map=c("genetic","physical"), placeUsing=c("qtl","correlation"), threshold=3, chr, use.orderMarkers=FALSE, verbose=FALSE, debugMode=0)
	
}

\arguments{
\item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details. }
\item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. If not supported, it will be created using data stored in population}
  \item{map}{ 
  Which map should be used for comparison:
  \itemize{
    \item{genetic}{ - genetic map from cross$maps$genetic.}
    \item{physical}{ - physical map from cross$maps$physical.}
  }
  }
\item{placeUsing}{ 
  How position of the new markers on the saturated map should be determinate:
  \itemize{
    \item{qtl}{ - placed between two markers with highest .}
    \item{correlation}{ - physical map from cross$maps$physical.}
  }
  }
 \item{threshold}{ Specifies threshold for selecting/rejecting markers (see \link{markerPlacementPlot}).}
 \item{chr}{ Specifies a subset of chromosomes analysis will be carried on. If none, all the chromosomes will be saturated.}
 \item{use.orderMarkers}{Should markers on the newly created map be ordered using\code{\link[qtl]{orderMarkers}} function.}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information.}
 \item{verbose}{ Be verbose.}
}

\value{
  An object of class \code{\link{population}}. See \code{\link{create.population}} for details.
}

\details{
This function saturates existing map (stored in the population object) with markers derived form gene expression data
(provided inside cross or population. Correlation matrix between those two is made and based on it, new markers sre
being placed on the map.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
	###
	yeastPopulation <- scan.qtls(yeastPopulation,verbose=TRUE,map="physical",step=2)
	cross <- cross.saturate(yeastPopulation,map="physical",verbose=TRUE,debugMode=2)
}

\seealso{
  \itemize{
    \item{\code{\link{reorganizeMarkersWithin}}}{ - Apply new ordering on the cross object usign ordering vector.}
    \item{\code{\link{assignChrToMarkers}}}{ - Create ordering vector from chromosome assignment vector.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or chromosome assignment vector.}
    \item{\code{\link{reduceChromosomesNumber}}}{ - Number of routines to reduce number of chromosomes of cross object.}
    \item{\code{\link{markerPlacementPlot}}}{ - Plot showing how many markers will be selected for map saturation with different thresholds.}
  }
}

\keyword{manip}
