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
	plotMarkerDistribution(population,marker,nrDistributions,logarithmic=FALSE)
}

\arguments{
 \item{population}{ population type object, must contain parental phenotypic data}
 \item{marker}{ number or name of marker to be printed }
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
	data(yeastPopulation)
	plotMarkerDistribution(yeastPopulation,2,2)
}

\seealso{
  \code{\link{plotChildrenExpression}}
  \code{\link{plotParentalExpression}}
  \code{\link{plotMapComparison}}
}

\keyword{manip}
