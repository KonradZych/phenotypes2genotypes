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
 \item{population}{ Population type object, must contain parental phenotypic data.}
 \item{markers}{ Numbers of markers to be plotted. }
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
	setwd(paste(.Library,"pheno2geno/data",sep="/"))
	population <- readFiles()
	plotChildrenExpression(population)
}

\seealso{
  \code{\link{plotChildrenExpression}}
  \code{\link{plotMarkerDistribution}}
  \code{\link{plotMapComparison}}
}

\keyword{manip}
