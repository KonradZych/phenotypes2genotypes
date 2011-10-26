\name{Assign chromosomes}
\alias{assignChromosomes}
\alias{sumMajorityRule}
\alias{majorityRule}
\alias{correlationRule}


\title{Assign chromosomes.}

\description{
  Creates new map out of data stored in population object and ordering it using old map.
}

\usage{
	assignChromosomes(population, n.chr, map=c("none","genetic","physical"), how = c(correlationRule,majorityRule,sumMajorityRule), reOrder=FALSE, verbose=FALSE, debugMode=0)
}

\arguments{
 \item{population}{ object of class population}
 \item{n.chr}{ number of chromosomes expected}
 \item{map}{ 
  Which map should be used for comparison:
  \itemize{
    \item{none}{ - returns non-ordered cross object}
    \item{genetic}{ - genetic map from cross$maps$genetic}
    \item{physical}{ - physical map from cross$maps$physical}
  }
  }
  \item{how}{ 
  Which function should be used for assigning chromosomes:
  \itemize{
    \item{correlationRule}{ - assigning chromosome from new map to old ones using mean corelation\
      between their markers}
    \item{majorityRule}{ - for every chromosome in cross for every marker checks the marker it is
      having highest correlation with. Checks on which chromosome this marker is placed in old map. For each of
      new chromosomee, old chromosome with most markers with high correlation is assigned.}
    \item{sumMajorityRule}{ - for every chromosome in cross for every marker checks the marker it is
        having highest correlation with. Checks on which chromosome this marker is placed in old map. For each of
        new chromosomes one or more of chromosomes from old map will be represented. Function sums correlations for
        each pair of those and for every new chromosomes assigns old chromosome with highest cumulative cor.}
  }
  }
 \item{reOrder}{ if FALSE (default) function returns vector with new ordering (unless map=none is selected, then cross is being returned) is TRUE,
 cross is reordered using created vector and then returned
 }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
}

\value{
  an object of R/qtl class cross or vector containing new ordering
}

\details{
  }

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	### simulating data
	population <- fakePopulation(type="riself")
	population <- findDiffExpressed(population)
	population <- toGenotypes(population,genotype="simulated",proportion=c(50,50),orderUsing="map_genetic",treshold=0.01)
	cross <- assignChromosomes(population,10,map="physical",verbose=TRUE)

}

\seealso{
  \code{\link{readFiles}}
  \code{\link{findDiffExpressed}}
  \code{\link{toGenotypes}}
}

\keyword{manip}

