\name{RPpval}
\alias{showRPpval}
\alias{plotRPpval}

\title{RP p-values visualisation}

\description{
  Printing out/plotting p-values calculated by the findDiffExpressed function.
}

\usage{
	showRPpval(population,markers=1:10)
	plotRPpval(population,markers=1:10,treshold=0.01)
}

\arguments{
 \item{population}{ an object of class \code{population}. See \code{\link{createPopulation}} for details.}
 \item{markers}{ numbers of markers to be printed }
 \item{treshold}{ value on which horizontal line will be plotted}
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
	showRPpval(population)
	plotRPpval(population)

}

\seealso{
	\itemize{
    \item \code{\link[RankProd]{RP}} - Add description
    \item \code{\link{findDiffExpressed}} - Finding differentially expressed genes.
  }
}

\keyword{manip}
