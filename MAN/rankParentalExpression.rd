\name{readParentalExpression}
\alias{readParentalExpression}

\title{Product Rank parental expression data.}

\description{
  Using Product Rank algorithm to select genes that are differentially expressed between parents.
}

\usage{
	rankParentalExpression(expressionParental,groupLabels=c(0,0,1,1),verbose=FALSE,debugMode=0,...)
}

\arguments{
 \item{expressionParental}{ Matrix rows - markers, cols - individuals. }
 \item{groupLabels}{ list cointaining 0 and 1 specifing two parental groups }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
 \item{...}{ additional arguments for RP function }
}

\value{
  RP function output object.
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
  \code{\link{readParentalExpression}}
  \code{\link{filterParentalExpression}}
}

\keyword{manip}
