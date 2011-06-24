\name{preprocessData}
\alias{preprocessData}
\alias{groupLabels}

\title{RankProduct analysis.}

\description{
  Using Rank Product analysis to select differentially expressed genes.
}

\usage{
	preprocessData(population,groupLabels=c(0,0,1,1),verbose=FALSE,debugMode=0,...)
}

\arguments{
 \item{population}{ population type object, must contain parental phenotypic data.}
 \item{groupLabels}{ Specify which column of parental data belongs to group 0 and which to group 1.}
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information }
 \item{...}{ Additional arguments passed to RP function. }
}

\value{
  Object of type RP (see ?RP for description), saved into population$founders$RP.
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
	\dontrun{
	setwd(paste(.Library,"pheno2geno/data",sep="/"))
	population <- readFiles()
	population <- preprocessData(population)
	}
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{toGenotypes}}
}

\keyword{manip}
