\name{selectChildrenExpression}
\alias{selectChildrenExpression}
\title{Marker selection.}

\description{
  Selects only markers that are differentialy expressed between parents.
}

\usage{
	selectChildrenExpression <- function(expressionChildren,expressionParental,verbose=FALSE,debugMode=0)
}

\arguments{
 \item{expressionChildren}{ matrix rows - markers, cols - individuals }
 \item{expressionParental}{ object containing parental expression data }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
}

\value{
  Matrix rows - genes with significant difference of expreession between PARENTS, cols - children individuals.
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
  \code{\link{correctChildrenExpression}}
  \code{\link{selectChildrenExpression}}
  \code{\link{readChildrenGenotypes}}
}

\keyword{manip}
