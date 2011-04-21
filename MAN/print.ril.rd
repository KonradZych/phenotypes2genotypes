\name{print.ril}
\alias{print.ril}


\title{RankProduct analysis.}

\description{
  Using Rank Product analysis to select differentially expressed genes.
}

\usage{
	print.ril(x)
}

\arguments{
 \item{x}{ Ril type object, must contain parental phenotypic data.}
}

\value{
  Object of type RP, saved into ril$parental$RP.
}

\details{
  Ril object description should go here.
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
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{toGenotypes}}
}

\keyword{manip}
