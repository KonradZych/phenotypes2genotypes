\name{preprocessData}
\alias{preprocessData}
\alias{ril}
\alias{groupLabels}

\title{RankProduct analysis.}

\description{
  Using Rank Product analysis to select differentially expressed genes.
}

\usage{
	preprocessData(ril,groupLabels=c(0,0,1,1),verbose=FALSE,debugMode=0,...)
}

\arguments{
 \item{ril}{ Ril type object, must contain parental phenotypic data.}
 \item{groupLabels}{ Specify which column of parental data belongs to group 0 and which to group 1.}
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information }
 \item{...}{ Additional arguments passed to RP function. }
}

\value{
  Object of type RP, saved into ril$parental$RP.
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
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{toGenotypes}}
}

\keyword{manip}
