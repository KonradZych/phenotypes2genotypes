\name{parentalRoutine}
\alias{parentalRoutine}
\alias{parentalFile}
\alias{groupLabels}
\alias{treshold}

\title{Reading parental expression data routine.}

\description{
  Read in parental data, process is using RankProd and gives nice output.
}

\usage{
	parentalRoutine(parentalFile="Gene_parental.txt",groupLabels=c(0,0,1,1),treshold=0.01,verbose=FALSE,debugMode=0,...)
}

\arguments{
 \item{parentalFile}{ file containing parental expression data }
 \item{groupLabels}{ list cointaining 0 and 1 specifing two parental groups }
 \item{treshold}{ treshold value - if rank is below, difference in expression values between genes are significant }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
 \item{...}{ additional arguments for RP function }
}

\value{
  Matrix nrow = nr of significant genes, ncol = 2, cointains mean expression values for each group.
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
	expressionParental <- readParentalExpression()
}

\seealso{
  \code{\link{readParentalExpression}}
  \code{\link{rankParentalExpression}}
  \code{\link{filterParentalExpression}}
}

\keyword{manip}
