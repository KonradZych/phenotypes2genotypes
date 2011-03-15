\name{bestClustering}
\alias{bestClustering}
\alias{groups}
\alias{iterations}
\alias{use}

\title{Iterational clustering using kmeans.}

\description{
  Perform clustering using kmeans in specified iteriations number.
}

\usage{
bestClustering(genotypeMatrix, groups, iterations, use="r", flip=0, verbose=FALSE, debugMode=0)
}

\arguments{
 \item{genotypeMatrix}{ Matrix with genotype values with: columns, individuals and on the rows, markers }
 \item{groups}{ number of groups data should be splitted to }
 \item{iterations}{ how many times clustering should be performed }
 \item{use}{ r -> use recombination matrix for clustering, c -> corelation matrix }
 \item{flip}{ specifies whether one of the rows that are being compared should be flipped(1) or not(0) }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
}

\value{
  A square matrix with numberes (0-iterations)
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
	library(qtl)
  	data(genotypes)
	data(expression_ratio)
	genotypesToCross(genotypes,expression_ratio)
}

\seealso{
  TODO
}

\keyword{manip}
