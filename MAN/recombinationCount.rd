\name{recombinationCount}
\alias{recombinationCount}

\title{Recombination factor count.}

\description{
  Counting number of recombination events between two markers.
}

\usage{
recombinationCount <- function(genotypicMatrix,flip=0,verbose=FALSE,debugMode=0)
}

\arguments{
 \item{genotypeMatrix}{ Matrix with genotype values with: columns, individuals and on the rows, markers }
 \item{flip}{ specifies whether one of the rows that are being compared should be flipped(1) or not(0) }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
}

\value{
  A square matrix with numberes (0-nr of individuals)
}

\details{
  TODO
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
	Under tender patronage of: Danny Arends \email{Danny.Arends@gmail.com}
}

\examples{
  	data(genotypes)
	reco <- recombinationCount(genotypes)
	#checks
	if(reco[1]!=0)	stop("Element not equal to expected\n")
	if(reco[34]!=30)	stop("Element not equal to expected\n")
	if(reco[445]!=24)	stop("Element not equal to expected\n")
	if(sum(dim(genotypes)) != 240)	stop("Wrong dimensions\n")
}

\seealso{
  TODO
}

\keyword{manip}
