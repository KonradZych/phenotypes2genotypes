\name{cross.saturate}
\alias{cross.saturate}

\title{Saving gff files.}

\description{
  Saving gff files.
}

\usage{
saveGff(cross, population, gffFile="population", verbose=FALSE)
}

\arguments{
\item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details. }
\item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. If not supplied, it will be created using data from the population object }
\item{gffFileCore}{Name of the gff files core where a physical location of the markers is stored for use in genome viewers. Four files will be saved - one with only
non-reduntant markers, one with redundat markers included and two with recombination breakpoints between the all/nonredundant markers.}
 \item{debugMode}{ Either use 1 or 2, this will modify the amount of information returned to the user. 1) Print out checks, 2) Print additional time information.}
 \item{verbose}{ Be verbose.}
}

\value{
  None.
}

\details{
This function saves gff files, that can be visualised using most of the genome viewers. The files contain physical location of markers and recombination breakpoints.
Therefore, physical map should be stored in an object of class population. Redundant markers are the markers having the same loaction on the genomic map, but different
on the physical map. These may be produced e.g. by cross.saturate function (markers that have QTL exactly on the position where an original marker is located). Also, as a
result of smoothing genotypping errors some markers may be put on the same position on the genetic map.
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
