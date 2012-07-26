\name{plotChildrenExpression}
\alias{plotChildrenExpression}
\alias{markers}

\title{Plotting routine for children expression data.}

\description{
  Plots offspring gene expression data in comparison with founders data.
}

\usage{
	plotChildrenExpression(population, markers=1:100)
}

\arguments{
 \item{population}{ An object of class \code{\link{population}}. See \code{\link{read.population}} for details. }
 \item{markers}{ Numbers of markers to be plotted.}
}

\value{
	None.
}

\details{
 Plots offspring expression data (boxplot) max value of parental expression (red triangle) min (blue triangle) and mean(line) for selected markers.
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
  \itemize{
    \item{\code{\link{plotParentalExpression}}}{ -  Plotting routine for parental gene expression data.}
    \item{\code{\link{plotMarkerDistribution}}}{ -  Plotting gene expression data for a single marker.}
}
}

\keyword{manip}
