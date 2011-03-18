\name{correctChildrenExpression}
\alias{correctChildrenExpression}


\title{Correction for children data.}

\description{
  Used if there are artificially iintroduced groups in data, e.g. data harvested in different conditions.
}

\usage{
	correctChildrenExpression <- function(expressionChildren,genotypeMatrix,verbose=FALSE,debugMode=0)
}

\arguments{
 \item{expressionChildren}{ matrix rows - markers, cols - individuals }
 \item{genotypeMatrix}{ matrix rows - markers, cols - individuals.}
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
  \code{\link{childrenRoutine}}
  \code{\link{selectChildrenExpression}}
  \code{\link{readChildrenGenotypes}}
}

\keyword{manip}
