\name{write.population}
\alias{write.population}

\title{Saving population object to file.}

\description{
  Saving object of class population to text file.
}

\usage{
	write.population(population,outputFile="population.txt",verbose=FALSE)
}

\arguments{
 \item{population}{ An object of class \code{population}. See \code{\link{createPopulation}} for details.}
 \item{outputFile}{ name of the output file}
 \item{verbose}{ be verbose}
}

\value{
  None.
}

\details{
  TODO
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	\dontrun{
		population <- fakePopulation()
		write.population(population,verbose=TRUE)
	}
}

\seealso{
	\itemize{
    \item \code{\link[RankProd]{RP}} - Add description
    \item \code{\link{readFiles}} - Add description
    \item \code{\link{findBiomarkers}} - Creating genotypes from children phenotypes
  }
}

\keyword{manip}
