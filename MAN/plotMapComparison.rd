\name{plotMapComparison}
\alias{plotMapComparison}
\alias{coloringMode}

\title{Plotting routine for comparison of two genetic maps.}

\description{
  Plotting routine for comparison of two genetic maps.
}

\usage{
	plotMapComparison(cross, coloringMode=1)
}

\arguments{
 \item{cross}{ R/qtl cross type object.}
 \item{coloringMode}{ 1 - rainbow colors 2 - black for cis and red for trans located markers. }
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
	ril <- readFiles()
	ril <- preprocessData(ril)
	#cross <- toGenotypes(ril,use="simulated",minChrLength=0,treshold=1,margin=50,max.rf=1,min.lod=0)
	#plotMapComparison(cross)
}

\seealso{
  \code{\link{plotChildrenExpression}}
  \code{\link{plotParentalExpression}}
}

\keyword{manip}
