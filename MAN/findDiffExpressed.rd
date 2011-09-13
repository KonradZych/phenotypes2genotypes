\name{findDiffExpressed}
\alias{findDiffExpressed}
\alias{groupLabels}

\title{Finding differentially expressed genes.}

\description{
  Using Rank Product analysis to select differentially expressed genes.
}

\usage{
	findDiffExpressed(population,verbose=FALSE,debugMode=0,...)
}

\arguments{
 \item{population}{ An object of class \code{population}. See \code{\link{createPopulation}} for details.}
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
	population <- findDiffExpressed(population)

}

\seealso{
	\itemize{
    \item \code{\link[RankProd]{RP}} - Add description
    \item \code{\link{readFiles}} - Add description
    \item \code{\link{toGenotypes}} - Creating genotypes from children phenotypes
    \item \code{\link{showRPpval}} - Printing out p-values calculated by the findDiffExpressed function
    \item \code{\link{plotRPpval}} - Plotting p-values calculated by the findDiffExpressed function
  }
}

\keyword{manip}
