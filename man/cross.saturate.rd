\name{cross.saturate}
\alias{cross.saturate}

\title{Saturate an existing genetic map.}

\description{
  Saturating an existing genetic map using markers derived from phenotype data.
}

\usage{
cross.saturate(population, cross, map=c("genetic","physical"), placeUsing=c("qtl","correlation"), flagged = c("remove","warn","ignore"), threshold=3,
  chr, env, use.orderMarkers=FALSE, verbose=FALSE, debugMode=0)

}

\arguments{
\item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details. }
\item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. If not supplied, it will be created using data from the population object }
  \item{map}{ 
  Which map should be used for comparison:
  \itemize{
    \item{genetic}{ - genetic map from cross$maps$genetic.}
    \item{physical}{ - physical map from cross$maps$physical.}
  }
  }
\item{placeUsing}{ 
  How should the position of the new markers on the saturated map be determinated:
  \itemize{
    \item{qtl}{ - position the new markers between / next to markers with high LOD score (see threshold). }
    \item{correlation}{ - position the new markers on the locations with the highest correction to markers on the physical map from cross$maps$physical. }
  }
  }
\item{flagged}{ 
  How to handle the markers influenced by epistatic or environmental interactions:
  \itemize{
    \item{remove}{ - warn about every marker affected and remove them. }
    \item{warn}{ - warn about every marker affected but leave them in. }
    \item{ignore}{ - leave them in. }
  }
  }
 \item{threshold}{ Specifies a threshold for the selection of new phenotype markers (see \link{markerPlacementPlot}).}
 \item{chr}{ When specified the algorithm only saturates a subset of chromosomes. If not specified, all the chromosomes will be saturated. }
 \item{env}{ Vector of environmental conditions - for each of the individuals specifies a condition. Ignored if missing.}
 \item{use.orderMarkers}{ If true the algorithm (after initial saturation) performs an \code{\link[qtl]{orderMarkers}} on the newly created map.}
 \item{verbose}{ Be verbose.}
 \item{debugMode}{ Either use 1 or 2, this will modify the amount of information returned to the user. 1) Print out checks, 2) Print additional time information.}
}

\value{
  An object of class \code{\link{population}}. See \code{\link{create.population}} for details.
}

\details{
This function saturates an existing map with markers derived from the phenotype data provided inside either the cross or population object. 
A correlation matrix between those two sets of markers is made, and new markers are assigned to the 'optimal' location on the map.
}

\author{
	Konrad Zych \email{k.zych@rug.nl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{k.zych@rug.nl}
}

\examples{
	data(yeastPopulation)
	cross <- cross.saturate(yeastPopulation,map="physical",verbose=TRUE,debugMode=2)
}

\seealso{
  \itemize{
    \item{\code{\link{reorganizeMarkersWithin}}}{ - Apply new ordering on the cross object usign ordering vector.}
    \item{\code{\link{assignChrToMarkers}}}{ - Create ordering vector from chromosome assignment vector.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or chromosome assignment vector.}
    \item{\code{\link{reduceChromosomesNumber}}}{ - Functions to reduce the number of chromosomes in a cross object.}
    \item{\code{\link{markerPlacementPlot}}}{ - Plot showing how many markers will be selected for map saturation with different thresholds.}
  }
}

\keyword{manip}
