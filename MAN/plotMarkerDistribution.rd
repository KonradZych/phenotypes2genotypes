\name{plotMarkerDistribution}
\alias{plotMarkerDistribution}
\alias{logarithmic}
\alias{phenotypeRow}

\title{plotMarkerDistribution}

\description{
  Plotting histogram out of gene expression data for a single marker and fitting specified number of normal distribution curves, using EM algorithm.
}

\usage{
	plotMarkerDistribution(population,marker,nrDistributions,logarithmic=FALSE)
}

\arguments{
 \item{population}{ an object of class \code{\link{population}}}
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
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
	plotMarkerDistribution(yeastPopulation,2,2)
}

\seealso{
  \code{\link{plotParentalExpression}} - plotting routine for parental gene expression data
  \code{\link{plotChildrenExpression}} - plotting routine for children gene expression data
}

\keyword{manip}
