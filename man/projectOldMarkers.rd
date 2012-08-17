\name{projectOldMarkers}
\alias{projectOldMarkers}

\title{Plotting routine which shows where markers from original map are located on saturated map.}

\description{
  Plotting routine which shows where markers from original map are located on saturated map.
}

\usage{
	projectOldMarkers(cross,population,map=c("genetic","physical"),label=c("positions","names","no"),...)
}

\arguments{
\item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
\item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details. }
 \item{map}{ Which map (from the ones stored in the population$maps) should be used to assigning chromosomes on the created map}
 \item{label}{ Should the old markers be labeled in the plot (options: position, name or off).}
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
