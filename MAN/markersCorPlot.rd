\name{markersCorPlot}
\alias{markersCorPlot}

\title{Plotting correlation between markers on two maps.}

\description{
  Plotting correlation between two maps together with markers placement (comparison of coverage).
}

\usage{
	markersCorPlot(cross, population, map=c("genetic","physical"), cmBetween=25, comparisonMethod = c(sumMajorityCorrelation,majorityCorrelation,meanCorrelation), chr, show.legend=FALSE, verbose=TRUE)
}

\arguments{
 \item{cross}{ R/qtl cross type object.}
  \item{population}{ An object of class \code{\link{population}}.}
 \item{map}{ Which map (from ones stored in population$maps) should be used fo assigning chromosomes on the created map.}
 \item{cmBetween}{ Offset between chromosomes specified in cM.}
 \item{comparisonMethod}{ Method used to compare chromosomes from the new map to the original ones while assigning:
   \itemize{
    \item{sumMajorityCorrelation}{ - For each chromosome in cross for every marker checks the marker it is
   having highest correlation with. Checks on which chromosome this marker is placed in old map. For each of
   new chromosomes one or more of chromosomes from old map will be represented. Function sums correlations for
   each pair of those and for every new chromosomes assigns old chromosome with highest cumulative cor.}
    \item{majorityCorrelation}{ - For each chromosome in cross for every marker checks the marker it is
   having highest correlation with. Checks on which chromosome this marker is placed in old map. For each of
   new chromosomee, old chromosome with most markers with high correlation is assigned.}
    \item{meanCorrelation}{ - Assigning chromosome from new map to old ones using sum of the mean correlation between their markers.}
    \item{majorityOfMarkers}{ - For each chromosome in the cross object (either created inside the function or provided
  by user) chromosome from original map, where most markers from new chromosome are is assigned.}
  }
 }
 \item{chr}{ Specifies subset of chromosomes to be shown.}
 \item{show.legend}{ Shall the legend be shown on the plot.}
 \item{verbose}{ Be verbose.}
}

\details{
Plots markers from moth old and new map as points and in the background - comparison between them done using selected comparison method.
}

\value{
	Plot.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
	data(yeastCross)
	markersCorPlot(yeastCross,yeastPopulation,map="physical")
}

\seealso{
  \itemize{
    \item{\code{\link{plotMapComparison}}}{ -  Plotting routine for comparison of two genetic maps.}
    \item{\code{\link{projectOldMarkers}}}{ -   Plotting routine for showing how markers from original map are placed on saturated map.}
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or chromosome assignment vector.}
}
}

\keyword{manip}
