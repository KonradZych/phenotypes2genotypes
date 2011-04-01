\name{toGenotypes}
\alias{toGenotypes}
\alias{use}
\alias{treshold}
\alias{overlapInd}
\alias{proportion}
\alias{margin}

\title{Creating genotypes from children phenotypes.}

\description{
  Creating genotypes from children phenotypes using parental data and saving cross object.
}

\usage{
	toGenotypes(ril, use=c("real","simulated","map"), treshold=0.01, overlapInd = 0, proportion = 50, margin = 15, verbose=FALSE, debugMode=0)
}

\arguments{
 \item{ril}{ Ril type object, must contain parental phenotypic data.}
 \item{treshold}{ If Rank Product pval for gene is lower that this value, we assume it is being diff. expressed.}
 \item{use}{ Which genotypic matrix should be saved to file, real - supported by user and read from file, simulated - made by toGenotypes, map - simulated data orderd using gff map}
 \item{overlapInd}{ Number of individuals that are allowed in the overlap }
 \item{proportion}{ Proportion of individuals expected to carrying a certain genotype }
 \item{margin}{ Proportion is allowed to varry between this margin (2 sided) }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
}

\value{
  Cross type object.
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
	setwd(paste(.Library,"pheno2geno/data",sep="/"))
	ril <- readFiles()
	ril <- preprocessData(ril)
	cross <- toGenotypes(ril,use="simulated")
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{preprocessData}}
}

\keyword{manip}
