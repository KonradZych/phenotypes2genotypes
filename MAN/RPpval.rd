\name{RPpval}
\alias{show.RPpval}
\alias{plot.RPpval}

\title{RP p-values visualisation}

\description{
  Printing out/plotting p-values calculated by the find.diff.expressed function.
}

\usage{
	show.RPpval(population,markers=1:10)
	plot.RPpval(population,thresholdRange=c(0.01,0.1,0.01))
}

\arguments{
\item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details. }
 \item{markers}{ Numbers of markers to be printed}
 \item{thresholdRange}{ Specifies in which range threshold will be checked (start, stop, step).}
}

\value{
  An object of class population, (see \code{\link{create.population}} for more details) with object of class \link[RankProd]{RP} saved into population$founders$RP.
}

\details{
  Those are two helper functions of \code{\link{find.diff.expressed}}. show.RPpval is printing to the screen p-values for specified markers, while plot.RPpval is showing
  how many markers will be selected using different thresholds.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{

	data(yeastPopulation)
	show.RPpval(yeastPopulation)
	plot.RPpval(yeastPopulation)

}

\seealso{
	\itemize{
    \item{\code{\link[RankProd]{RP}}}{ - Perform rank product method to identify differentially expressed genes.}
    \item{\code{\link{find.diff.expressed}}}{ - Select differentially expressed genes using Rank Product or student t-test analysis.}
    \item{\code{\link{generateBiomarkers}}}{ - Creating genotypes from children phenotypes.}
    \item{\code{\link{show.RPpval}}}{- Printing out p-values calculated by the find.diff.expressed function.}
    \item{\code{\link{plot.RPpval}}}{ - Plotting p-values calculated by the find.diff.expressed function.}
  }
}

\keyword{manip}
