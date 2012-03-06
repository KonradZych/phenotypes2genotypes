\name{projectOldMarkers}
\alias{projectOldMarkers}

\title{Projecting old markers on saturated map.}

\description{
  Plotting routine for showing how markers from original map are placed on saturated map.
}

\usage{
	projectOldMarkers(cross,population,map=c("genetic","physical"),label=c("positions","names","no"),...)
}

\arguments{
\item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
\item{population}{ An object of class \code{\link{population}}. See \code{\link{createPopulation}} for details. }
 \item{map}{ Which map (from ones stored in population$maps) should be used fo assigning chromosomes on the created map}
 \item{label}{ How old markers should be labeled on the plot, with their positions, their names or none.}
 \item{...}{Parameters passed to \code{\link[qtl]{plot.qtl}}.}
}

\value{
	None.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
	data(yeastCross)
	projectOldMarkers(yeastCross,yeastPopulation,map="physical")
}

\seealso{
  \itemize{
    \item{\code{\link{plotMapComparison}}}{ -  Plotting routine for comparison of two genetic maps.}
    \item{\code{\link{markersCorPlot}}}{ -  Plotting correlation between two maps together with markers placement (comparison of coverage).}
}
}

\keyword{manip}
