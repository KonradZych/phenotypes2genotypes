\name{plotParentalExpression}
\alias{plotParentalExpression}

\title{Plotting routine for parental expression data.}

\description{
  Plots parental data in two colors (two parental groups) and mean of values for each marker.
}

\usage{
	plotParentalExpression(population, markers=1:100, groupLabels=c(0,0,1,1))
}

\arguments{
 \item{population}{ population type object, must contain parental phenotypic data.}
 \item{markers}{ Numbers of markers to be plotted. }
 \item{groupLabels}{ Specify which column of parental data belongs to group 0 and which to group 1.}
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
	### simulating data
	population <- fakePopulation()
	### plotting 
	plotParentalExpression(population)
}

\seealso{
  \code{\link{plotChildrenExpression}}
  \code{\link{plotMarkerDistribution}}
  \code{\link{plotMapComparison}}
}

\keyword{manip}
