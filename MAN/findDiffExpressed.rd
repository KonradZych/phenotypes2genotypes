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
 \item{population}{ An object of class \code{population}. See \code{\link{createPopulation}} for details.}
 \item{groupLabels}{ Specify which column of parental data belongs to group 0 and which to group 1.}
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information }
 \item{...}{ Additional arguments passed to RP function. }
}

\value{
  Object of class population, (see \code{\link{createPopulation}} for more details about structure) with object of class \link[RankProd]{RP} saved into population$founders$RP.
}

\details{
  TODO
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{

	population <- fakePopulation()
	population <- preprocessData(population,c(0,0,0,1,1,1))

}

\seealso{
	\itemize{
    \item \code{\link[RankProd]{RP}} - Add description
    \item \code{\link{readFiles}} - Add description
    \item \code{\link{toGenotypes}} - Creating genotypes from children phenotypes
  }
}

\keyword{manip}
