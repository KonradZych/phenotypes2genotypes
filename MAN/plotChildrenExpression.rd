\name{plotChildrenExpression}
\alias{plotChildrenExpression}
\alias{markers}

\title{Plotting routine for children expression data.}

\description{
  Plots children expression data (boxplot) max value of parental exptression (red trangle) min (blue triangle) and mean(line) for selected markers.
}

\usage{
	plotChildrenExpression(population, markers=1:100)
}

\arguments{
  \item{population}{ an object of class \code{\link{population}}}
 \item{markers}{ Numbers of markers to be plotted. }
}

\value{
	Plot.
}

\details{
  TODO
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
	### plotting only 10 markers for clearer image
	plotChildrenExpression(yeastPopulation,10:20)

}

\seealso{
  \code{\link{plotParentalExpression}} - plotting routine for parental gene expression data
  \code{\link{plotMarkerDistribution}} - plotting gene expression data for a single marker
}

\keyword{manip}
