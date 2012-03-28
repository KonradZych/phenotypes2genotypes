\name{find.mixups}
\alias{find.mixups}

\title{Find sample mix-ups}

\description{
  Finding possible sample mix-ups in the data.
}

\usage{
  find.mixups(population,map=c("genetic","physical"),n.qtls=50,threshold=15,verbose=FALSE)
}

\arguments{
  \item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details.}
   \item{map}{ Which map ( from ones stored in population$maps) contains information about positions of original markers.}
 \item{n.qtls}{ Number of qtls that should be scanned.}
 \item{threshold}{ How big percentage of being flagged must be for the individual to be considered as a mix-up.}
 \item{verbose}{ Be verbose.}
}

\value{
  An vactor with scores for each of the individuals. 
}

\details{
	This function scans data, looking for qtls. When requested number of qtls are found, they are being checked for outliers. If a value for an
  individual is not in the range (mean - 3*sd, mean+3*sd) of its group, it's being flagged. For each of the markers, percentage of being flagged is returned.
  Additionally, individuals crossing threshold set by the user are being printed out with the warning about possible mix-up.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
  scores <- find.mixups(yeastPopulation,map="physical",n.qtls=10,threshold=5,verbose=FALSE)
  plot(scores[[2]])
}

\seealso{
  \itemize{
    \item{\code{\link{read.population}}}{ - Load genotype, phenotype, genetic map data files into R environment into a population object.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or vector showing how chromosomes should be assigned.}
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{find.diff.expressed}}}{ - Using Rank Product or student t-test analysis to select differentially expressed genes.}
  }
}

\keyword{manip}
