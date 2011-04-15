\name{cleanMap}
\alias{cleanMap}


\title{Creating genotypes from children phenotypes.}

\description{
  Creating genotypes from children phenotypes using parental data and saving cross object.
}

\usage{
	cleanMap(cross, difPercentage, verbose=FALSE, debugMode=0)
}

\arguments{
 \item{cross}{ R/qtl cross type object.}
 \item{difPercentage}{ If removing marker will make map shorter by this percentage of its length, then it will be dropped.}
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
	#cross <- toGenotypes(ril,use="simulated",minChrLength=0,treshold=0.5,margin=50,max.rf=10)
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{preprocessData}}
  \code{\link{toGenotypes}}
}

\keyword{manip}
