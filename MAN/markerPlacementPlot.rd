\name{markerPlacementPlot}
\alias{markerPlacementPlot}


\title{Plot number of markers selected.}

\description{
 Plot number of markers selected with different thresholds.
}

\usage{
	markerPlacementPlot(population, placeUsing=c("qtl","correlation"), cross)
}

\arguments{
\item{population}{ An object of class \code{\link{population}}. See \code{\link{createPopulation}} for details. }
\item{placeUsing}{ 
  How position of the new markers on the saturated map should be determinate:
  \itemize{
    \item{qtl}{ - placed between two markers with highest .}
    \item{correlation}{ - physical map from cross$maps$physical.}
  }
  }
\item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
}

\value{
  None.
}

\details{
This plot is really useful while saturating existing map (using \code{\link{cross.saturate}}). It helps choose best threshold for marker selection, showing how much markers will
be selected with different threshold values.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastCross)
	data(yeastPopulation)
	markerPlacementPlot(yeastPopulation,placeUsing="correlation",yeastCross)
	markerPlacementPlot(yeastPopulation,placeUsing="qtl",yeastCross)
}

\seealso{
  \itemize{
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{reorganizeMarkersWithin}}}{ - Apply new ordering on the cross object usign ordering vector.}
    \item{\code{\link{assignedChrToMarkers}}}{ - Create ordering vector from chromosome assignment vector.}
    \item{\code{\link{reduceChromosomesNumber}}}{ - Number of routines to reduce number of chromosomes of cross object.}
    \item{\code{\link{findBiomarkers}}}{ - Creating genotype markers out of gene expression data.}
}
}

\keyword{manip}
