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
  \item{map}{ Which map should be used to determine the ordering / positions of original markers.}
  \item{n.qtls}{ Number of qtls that we use for scanning for mix-ups. }
  \item{threshold}{ When an individual is not matching the expected genotype more the x \% of the time. The individual should be considered as being a mix-up. }
  \item{verbose}{ Be verbose.}
}

\value{
  An vector with for each individual a percentage that shows how many times an individual didn't match the expected genotype.
}

\details{
	After scanning for the requested number of QTLs, each individual is checked if their genotype is matching the expected genotype. If an 
  individuals expression value is not in the range of the expected genotype (mean - 3*sd, mean+3*sd), it's receives a penaltie.
  After which the individuals above the threshold are being printed with a warning about possible mix-up.
}

\author{
	Konrad Zych \email{k.zych@rug.nl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{k.zych@rug.nl}
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
