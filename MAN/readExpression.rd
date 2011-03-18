\name{readParentalExpression}
\alias{readParentalExpression}

\title{Reading parental expression data from file.}

\description{
  Read in parental data.
}

\usage{
	readParentalExpression(parentalFile="Gene_parental.txt",verbose=FALSE,debugMode=0)
}

\arguments{
 \item{parentalFile}{ file containing parental expression data }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
}

\value{
  Matrix rows - markers, cols - individuals.
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
	
}

\seealso{
  \code{\link{parentalRoutine}}
  \code{\link{rankParentalExpression}}
  \code{\link{filterParentalExpression}}
}

\keyword{manip}
