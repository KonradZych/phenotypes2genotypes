\name{read.population}
\alias{read.population}

\title{Readinf population object from file.}

\description{
  Loading object of class population from text file.
}

\usage{
	read.population(filename="population.txt",verbose=FALSE)
}

\arguments{
 \item{filename}{ name of the population file}
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
	read.population("population.txt")
	}
}

\seealso{
	\itemize{
    \item \code{\link[RankProd]{RP}} - Add description
    \item \code{\link{readFiles}} - Add description
    \item \code{\link{toGenotypes}} - Creating genotypes from children phenotypes
  }
}

\keyword{manip}
