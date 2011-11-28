\name{Create new map}
\alias{createNewMap}
\alias{assignMaximumNoConflicts}
\alias{assignMaximum}


\title{Creating de novo genetic map.}

\description{
  Creating de novo genetic map.
}

\usage{
	createNewMap(population, cross, n.chr, map=c("none","genetic","physical"), comparisonMethod = c(sumMajorityCorrelation,majorityCorrelation,meanCorrelation,majorityOfMarkers), 
assignFunction=c(assignMaximumNoConflicts,assignMaximum), reOrder=TRUE, use.orderMarkers=FALSE, verbose=FALSE, debugMode=0)
	
}

\arguments{
 \item{population}{ an object of class \code{\link{population}}}
 \item{cross}{ an object of R/qtl class cross if missing, it will be created from population object}
 \item{n.chr}{number of chromosomes expected on the map}
 \item{map}{ which map (from ones stored in population$maps) should be used fo assigning chromosomes on the created map}
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
 \item{use.orderMarkers}{should markers on the newly created map be ordered using R/qtl orderMarkers funtion}
 \item{verbose}{ be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
}

\value{
  an object of class cross or vector showing how chromosomes should be assigned, to be used with \code{\link{assignedChrToMarkers}}
}

\details{
createNewMap function creates new genetic map using genotypes simulated by findBiomarkers function. Then it uses information provided by user to
assign number to newly created chromosomes.
}
\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
	cross <- createNewMap(yeastPopulation,n.chr=16,verbose=TRUE,map="physical",comparisonMethod=sumMajorityCorrelation, use.orderMarkers=FALSE)
}

\seealso{
  \code{\link{reorganizeMarkersWithin}} - apply new ordering on the cross object usign ordering vector
  \code{\link{assignedChrToMarkers}} - create ordering vector from chromosome assignment vector
  \code{\link{saturateExistingMap}} - saturate existing map
  \code{\link{reduceChromosomesNumber}} - number of routines to reduce number of chromosomes of cross object
  \code{\link{findBiomarkers}} - creating genotype markers  out of gene expression data.
}

\keyword{manip}
