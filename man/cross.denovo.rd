\name{cross.denovo}
\alias{cross.denovo}

\title{Create a de novo genetic map from a population object.}

\description{
  Create a de novo genetic map from offspring phenotype data stored in a population object
}

\usage{
cross.denovo(population, 
             n.chr, 
             orderUsingMap=FALSE, 
             map=c("none", "genetic", "physical"), 
             comparisonMethod = c(sumMajorityCorrelation, majorityCorrelation, meanCorrelation, majorityOfMarkers),
             assignFunction=c(assignMaximumNoConflicts, assignMaximum), 
             reOrder=TRUE, 
             use.orderMarkers=FALSE, 
             cross, 
             verbose=FALSE, 
             debugMode=0)
	
}

\arguments{
\item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details. }
 \item{n.chr}{ Number of chromosomes expected on the map.}
 \item{orderUsingMap}{ Shall markers in the result cross be ordered using one of the maps in population object.}
 \item{map}{ Which map ( from ones stored in population$maps) should be used for assigning chromosomes on the created map.}
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
 \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
 \item{verbose}{ be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
}

\value{
  When reordering this will produce an object of class \code{cross}, otherwise (reOrder=FALSE) 
  a chromosomes assignment vector (See \code{\link{assignChrToMarkers}} ) is produced which can be used to manual reorder the markers.
}

\details{
cross.denovo function creates new genetic map using genotypes simulated by \code{\link{generate.biomarkers}} function. Then it uses information provided by user to
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
    \item{\code{\link{assignChrToMarkers}}}{ - Create ordering vector from chromosome assignment vector.}
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{reduceChromosomesNumber}}}{ - Number of routines to reduce number of chromosomes of cross object.}
    \item{\code{\link{generate.biomarkers}}}{ - Creating genotype markers out of gene expression data.}
}
}

\keyword{manip}
