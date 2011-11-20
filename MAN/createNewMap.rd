\name{Create new map}
\alias{createNewMap}
\alias{assignMaximumNoConflicts}
\alias{assignMaximum}


\title{Creating de novo genetic map.}

\description{
  Creating de novo genetic map.
}

\usage{
	createNewMap(population, n.chr, map=c("none","genetic","physical"), comparisonMethod = c(sumMajorityCorrelation,majorityCorrelation,meanCorrelation,majorityOfMarkers), 
	assignFunction=c(assignMaximumNoConflicts,assignMaximum),reOrder=TRUE, use.orderMarkers=TRUE, verbose=FALSE, debugMode=0)
	
}

\arguments{
 \item{population}{ an object of class population}
 \item{n.chr}{number of chromosomes expected on the map}
 \item{map}{ which map (from ones stored in population$maps) should be used fo assigning chromosomes on the created map
 \item{comparisonMethod}{method used tocompare chromosomes from the new map to the original ones while assigning
   \itemize{
    \item{\code{\link{sumMajorityCorrelation}}}{}
    \item{\code{\link{majorityCorrelation}}}{}
    \item{\code{\link{meanCorrelation}}}{}
    \item{\code{\link{majorityOfMarkers}}}{}
  }
 }
 \item{assignFunction}{function used to assign chromosomes on the created map, in both cases for every chromosome from the new map, original chromosome with maximal score is assigned, but 
	if one of the original chromosomes is assigned to more then one of new ones:
	\itemize{
		\item{assignMaximumNoConflicts}{additional step is performed to make sure each of the original chromosomes is used only once}
		\item{assignMaximum}{those two are being merged}
	}
 }
 \item{reOrder}{ if TRUE, cross object is returned, FALSE - vector showing how chromosomes should be assigned}
 \item{use.orderMarker}{should markers on the newly created map be ordered using R/qtl orderMarkers funtion}
 \item{verbose}{ be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
}

\value{
  an object of class cross or vector showing how chromosomes should be assigned, to be used with \code{\link{reorganizeMarkersWithin}}
}

\details{
postProc function is makign use of most basic piece of inromation possible which is number of linkage groups
(normally equall to number of chromosome) expected. It is using \code{\link{formLinkageGroups}} functions from
R/qtl package with number of different parameter values and afterwards assesses which combination was the best one.
}
\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
	cross <- createNewMap(yeastPopulation,16,verbose=TRUE,map="physical",comparisonMethod=sumMajorityCorrelation, use.orderMarkers=FALSE)
}

\seealso{
  \code{\link{orderChromosomes}} - ordering chromosomes of an object of class cross using majority rule
  \code{\link{rearrangeMarkers}} - rearrangeing markers inside an object of class cross using correlation
  \code{\link{reduceChromosomesNumber}} - removing all but certain number of chromosomes from an object of class cross
}

\keyword{manip}
