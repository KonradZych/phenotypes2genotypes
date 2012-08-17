\name{write.population}
\alias{write.population}

\title{Writes a population object to file.}

\description{
  Writes a population object to file, for easy loading of intermediate data later.
}

\usage{
	write.population(population,outputFile="population.txt",verbose=FALSE)
}

\arguments{
 \item{population}{ An object of class \code{population}. See \code{\link{create.population}} for details. }
 \item{outputFile}{ Name of the output file. }
 \item{verbose}{ Be verbose. }
}

\value{
  None.
}

\details{
  This function writes an object of class \code{population} into a file.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	\dontrun{
		population <- fake.population()
		write.population(population,verbose=TRUE)
	}
}

\seealso{
  \itemize{
    \item{\code{\link{add.to.population}}}{ - Adding data to existing population object.}
    \item{\code{\link{create.population}}}{ - Create new object of class population.}
    \item{\code{\link{read.population}}}{ - Create new object of class population.}
  }
}

\keyword{manip}
