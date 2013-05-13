\name{markerPlacementPlot}
\alias{markerPlacementPlot}


\title{Plot number of markers selected.}

\description{
 Plot number of markers selected with different thresholds.
}

\usage{
  markerPlacementPlot(population, placeUsing=c("qtl","correlation"), 
    thrRange=c(1,5,1),cross,verbose=FALSE)
}

\arguments{
\item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details. }
\item{placeUsing}{ 
  How position of the new markers on the saturated map should be determinate:
  \itemize{
    \item{qtl}{ - placed between two markers with highest .}
    \item{correlation}{ - physical map from cross$maps$physical.}
  }
  }
\item{thrRange}{ Range of the threshold to be checked. Specified in a format(start,stop,step).}
\item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
 \item{verbose}{ Be verbose.}
}

\value{
  None.
}

\details{
This plot is really useful while saturating existing map (using \code{\link{cross.saturate}}). It helps choose best threshold for marker selection, showing how much markers will
be selected with different threshold values.
}

\author{
	Konrad Zych \email{k.zych@rug.nl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{k.zych@rug.nl}
}

\examples{
	data(yeastCross)
	data(yeastPopulation)
	markerPlacementPlot(yeastPopulation,placeUsing="qtl",cross=yeastCross)
	markerPlacementPlot(yeastPopulation,placeUsing="correlation",cross=yeastCross)
}

\seealso{
  \itemize{
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{reorganizeMarkersWithin}}}{ - Apply new ordering on the cross object usign ordering vector.}
    \item{\code{\link{assignChrToMarkers}}}{ - Create ordering vector from chromosome assignment vector.}
    \item{\code{\link{reduceChromosomesNumber}}}{ - Number of routines to reduce number of chromosomes of cross object.}
    \item{\code{\link{generate.biomarkers}}}{ - Creating genotype markers out of gene expression data.}
}
}

\keyword{manip}
