\name{plotMarkerDistribution}
\alias{plotMarkerDistribution}
\alias{logarithmic}
\alias{phenotypeRow}

\title{plotMarkerDistribution}

\description{
  Plotting histogram of distribution of values for single marker and specified number
 of normal distribution curves, fitted to data using EM algorithm
}

\usage{
	plotMarkerDistribution(phenotypeRow,nrDistributions,logarithmic=FALSE)
}

\arguments{
 \item{phenotypeRow}{ phenotypic data for single marker}
 \item{nrDistributions}{ numbers of normal distributions to be fitted}
 \item{logarithmic}{ TRUE - log(data) is used instead of raw data }
}

\value{
	Plot.
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
	#TODO
}

\seealso{
  \code{\link{plotChildrenExpression}}
  \code{\link{plotParentalExpression}}
  \code{\link{plotMapComparison}}
}

\keyword{manip}
