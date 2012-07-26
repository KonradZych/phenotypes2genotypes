\name{plotParentalExpression}
\alias{plotParentalExpression}

\title{Plotting routine for parental expression data.}

\description{
  Plots parental gene expression data.
}

\usage{
	plotParentalExpression(population, markers=1:100, groupLabels=c(0,0,1,1))
}

\arguments{
\item{population}{ An object of class \code{\link{population}}. See \code{\link{create.population}} for details. }
 \item{markers}{ Numbers of markers to be plotted. }
 \item{groupLabels}{ Specify which column of parental data belongs to group 0 and which to group 1.}
}

\value{
	None.
}

\details{
  Plots parental gene expression data in two colors (two parental groups) and mean of values for each marker.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	data(yeastPopulation)
	### plotting 
	plotParentalExpression(yeastPopulation)
}

\seealso{
  \itemize{
    \item{\code{\link{plotChildrenExpression}}}{ -  Plotting routine for children gene expression data.}
    \item{\code{\link{plotMarkerDistribution}}}{ -  Plotting gene expression data for a single marker.}
}
}

\keyword{manip}
