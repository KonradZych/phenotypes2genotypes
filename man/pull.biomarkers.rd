\name{pull.biomarkers}
\alias{pull.biomarkers}

\title{ Extract the detected biomarkers from a population object.}

\description{
  Extract the detected biomarkers from a population object, or select biomarkers that best match a certain pattern.
}

\usage{
  pull.biomarkers(population, pattern, verbose=FALSE)
}

\arguments{
  \item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details. }
  \item{pattern}{ Vector containg pattern to be matched in markers.}
  \item{verbose}{ Be verbose.}
}

\value{
  Matrix of all markers / vector with markers best matching the specified pattern.
}

\details{
	After running \code{\link{generate.biomarkers}} function, biomarkers are stored inside \code{\link{population}} class object.
  Use the pull.biomarkers function to extract them from the popuklation object into a matrix. This will return a matrix will 
  all the markers or when pattern is specified a vector with biomarkers best matching the pattern.
}

\author{
	Konrad Zych \email{k.zych@rug.nl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{k.zych@rug.nl}
}

\examples{
  data(yeastPopulation)
  markers <- pull.biomarkers(yeastPopulation,verbose=TRUE)
  bestMarker <- pull.biomarkers(yeastPopulation,round(runif(109)),verbose=TRUE)
}

\seealso{
  \itemize{
    \item{\code{\link{generate.biomarkers}}}{ - Create genotype markers  out of gene expression data.}
    \item{\code{\link{read.population}}}{ - Load genotype, phenotype, genetic map data files into R environment into a population object.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or vector showing how chromosomes should be assigned.}
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{find.diff.expressed}}}{ - Using Rank Product or student t-test analysis to select differentially expressed genes.}
  }
}

\keyword{manip}
