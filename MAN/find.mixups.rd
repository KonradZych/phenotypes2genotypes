\name{find.mixups}
\alias{find.mixups}

\title{Find sample mix-ups}

\description{
  Finding possible sample mix-ups in the data.
}

\usage{
  find.mixups(population,n.qtls=50,threshold=5,verbose=FALSE)
}

\arguments{
  \item{population}{ An object of class \code{\link{population}}. See \code{\link{createPopulation}} for details. }
 \item{n.qtls}{ Number of qtls that should be scanned.}
 \item{threshold}{ How many time and individual must be flagged to be considered a mix-up.}
 \item{verbose}{ Be verbose.}
}

\value{
  An vactor with scores for each of the individuals. 
}

\details{
	This function scans data, looking for qtls. When requested number of qtls are found, they are bieng checked for outliers. If a value for an
  individual is not in the range (mean - 3*sd, mean+3*sd) of its group, it's being flagged. Sum of flags for each of the individuals is returned
  in a for of a matrix. Additionally, individuals crossing threshold set by the user are being printed out with the warning about possible mix-up.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
  scores <- find.mixups(yeastPopulation,10)
  plot(scores)
}

\seealso{
  \itemize{
    \item{\code{\link{readFiles}}}{ - Load genotype, phenotype, genetic map data files into R environment into a population object.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or vector showing how chromosomes should be assigned.}
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{findDiffExpressed}}}{ - Using Rank Product or student t-test analysis to select differentially expressed genes.}
  }
}

\keyword{manip}
