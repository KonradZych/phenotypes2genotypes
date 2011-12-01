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
\item{population}{ An object of class \code{\link{population}}. See \code{\link{createPopulation}} for details. }
 \item{markers}{ Numbers of markers to be printed}
 \item{treshold}{ Value on which horizontal line will be plotted.}
}

\value{
  An object of class population, (see \code{\link{createPopulation}} for more details) with object of class \link[RankProd]{RP} saved into population$founders$RP.
}

\details{
  Those are two helper functions of \code{\link{findDiffExpressed}}. One is printing out, while the other is plotting, results of the analysis.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{

	data(yeastPopulation)
	showRPpval(yeastPopulation)
	plotRPpval(yeastPopulation)

}

\seealso{
	\itemize{
    \item{\code{\link[RankProd]{RP}}}{ - Perform rank product method to identify differentially expressed genes.}
    \item{\code{\link{findDiffExpressed}}}{ - Select differentially expressed genes using Rank Product or student t-test analysis.}
    \item{\code{\link{findBiomarkers}}}{ - Creating genotypes from children phenotypes.}
    \item{\code{\link{showRPpval}}}{- Printing out p-values calculated by the findDiffExpressed function.}
    \item{\code{\link{plotRPpval}}}{ - Plotting p-values calculated by the findDiffExpressed function.}
  }
}

\keyword{manip}
