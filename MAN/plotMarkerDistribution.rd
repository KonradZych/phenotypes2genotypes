\name{plotMarkerDistribution}
\alias{plotMarkerDistribution}
\alias{logarithmic}
\alias{phenotypeRow}

\title{plotMarkerDistribution}

\description{
  Plotting distribution of gene expression values of a single marker.
}

\usage{
	plotMarkerDistribution(population,marker,nrDistributions,logarithmic=FALSE)
}

\arguments{
\item{population}{ An object of class \code{\link{population}}. See \code{\link{createPopulation}} for details. }
 \item{marker}{ Number or name of the marker to be printed.}
 \item{nrDistributions}{ Number of normal distributions to be fitted.}
 \item{logarithmic}{ TRUE - log(data) is used instead of raw data.}
}

\value{
	None.
}

\details{
  Plotting histogram out of gene expression data for a single marker and fitting specified number of normal distribution curves, using EM algorithm.
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
  \itemize{
    \item{\code{\link{plotParentalExpression}}}{ -  Plotting routine for parental gene expression data.}
    \item{\code{\link{plotChildrenExpression}}}{ -   Plotting routine for children gene expression data.}
}
}

\keyword{manip}
