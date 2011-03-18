\name{filterParentalExpression}
\alias{filterParentalExpression}

\title{Filtering parental expression data.}

\description{
  Filters parental expression data, leaving only matrkers that show significant differencies in expression between parental groups.
}

\usage{
	filterParentalExpression(expressionParental,rankParental,groupLabels,treshold=0.01,verbose=FALSE,debugMode=0)
}

\arguments{
 \item{expressionParental}{ matrix rows - markers, cols - individuals }
 \item{rankParental}{ RP function output object }
 \item{groupLabels}{ list cointaining 0 and 1 specifing two parental groups }
 \item{treshold}{ treshold value - if rank is below, difference in expression values between genes are significant }
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
  \code{\link{readParentalExpression}}
}

\keyword{manip}
