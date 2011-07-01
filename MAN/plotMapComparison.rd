\name{plotMapComparison}
\alias{plotMapComparison}
\alias{coloringMode}

\title{Plotting routine for comparison of two genetic maps.}

\description{
  Plotting routine for comparison of two genetic maps.
}

\usage{
	plotMapComparison(cross, coloringMode=1, map=c("genetic","physical"))
}

\arguments{
 \item{cross}{ R/qtl cross type object.}
 \item{coloringMode}{ 1 - rainbow colors 2 - black for cis and red for trans located markers. }
 \item{map}{ 
  Which genotypic matrix should be saved to file:
  \itemize{
    \item{genetic}{ - genetic map from cross$maps$genetic}
    \item{physical}{ - physical map from cross$maps$physical}
  }
  }
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
	population <- fakePopulation()
	cross <- toGenotypes(population,genotype="real",orderUsing="map_genetic")
	plotMapComparison(cross,map="genetic")
}

\seealso{
  \code{\link{plotChildrenExpression}}
  \code{\link{plotParentalExpression}} 
}

\keyword{manip}
