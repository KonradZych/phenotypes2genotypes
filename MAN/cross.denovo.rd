\name{cross.denovo}
\alias{cross.denovo}
\alias{assignMaximumNoConflicts}
\alias{assignMaximum}
\alias{sumMajorityCorrelation}
\alias{majorityCorrelation}
\alias{meanCorrelation}
\alias{majorityOfMarkers}


\title{Creating de novo genetic map.}

\description{
  Creating de novo genetic map form gene expression data.
}

\usage{
cross.denovo(population, cross, n.chr, map=c("none","genetic","physical"), comparisonMethod = c(sumMajorityCorrelation,majorityCorrelation,meanCorrelation,majorityOfMarkers), 
assignFunction=c(assignMaximumNoConflicts,assignMaximum), reOrder=TRUE, use.orderMarkers=FALSE, verbose=FALSE, debugMode=0)
	
}

\arguments{
\item{population}{ An object of class \code{\link{population}}. See \code{\link{createPopulation}} for details. }
\item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
 \item{n.chr}{ Number of chromosomes expected on the map.}
 \item{map}{ Which map ( from ones stored in population$maps) should be used fo assigning chromosomes on the created map.}
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
  An object of class \code{cross} or vector showing how chromosomes should be assigned, to be used with \code{\link{assignedChrToMarkers}}
}

\details{
cross.denovo function creates new genetic map using genotypes simulated by generateBiomarkers function. Then it uses information provided by user to
assign number to newly created chromosomes.
}
\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
	cross <- cross.denovo(yeastPopulation,n.chr=16,verbose=TRUE,map="physical",comparisonMethod=sumMajorityCorrelation, use.orderMarkers=FALSE)
}

\seealso{
  \itemize{
    \item{\code{\link{reorganizeMarkersWithin}}}{ - Apply new ordering on the cross object usign ordering vector.}
    \item{\code{\link{assignedChrToMarkers}}}{ - Create ordering vector from chromosome assignment vector.}
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{reduceChromosomesNumber}}}{ - Number of routines to reduce number of chromosomes of cross object.}
    \item{\code{\link{generateBiomarkers}}}{ - Creating genotype markers out of gene expression data.}
}
}

\keyword{manip}
