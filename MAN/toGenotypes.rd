\name{toGenotypes}
\alias{toGenotypes}

\title{Convert a phenotypematrix into suitable genotypes}

\description{
  Convert a phenotypematrix into suitable genotypes
}

\usage{
toGenotypes(expressionMatrix, splitFUN = zero, overlapInd = 0, proportion = 50, margin = 5, genotypes = c(0,1), verbose=FALSE, debugMode=0)
}

\arguments{
 \item{expressionMatrix}{ Matrix with expression values with: columns, individuals and opn the rows, markers }
 \item{splitFUN}{ function used to split values }
 \item{overlapInd}{ Number of individuals that are allowed in the overlap }
 \item{proportion}{ Proportion of individuals expected to carrying a certain genotype }
 \item{margin}{ Proportion is allowed to varry between this margin (2 sided) }
 \item{genotypes}{ User defined genotypes for the output matrix }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
}

\value{
  A matrix with genotypes
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
  #TODO
}

\seealso{
  TODO
}

\keyword{manip}
