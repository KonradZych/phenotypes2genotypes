\name{childrenRoutine}
\alias{childrenRoutine}
\alias{parentalFile}
\alias{groupLabels}
\alias{treshold}

\title{Reading children expression data routine.}

\description{
  Read in children data, process is and gives nice output.
}

\usage{
	childrenRoutine(childrenFile="Gene_quant.txt",genotypeFile="Genotypes.txt",correction=TRUE,expressionParental=NULL,verbose=FALSE,debugMode=0)
}

\arguments{
 \item{childrenFile}{ file containing parental expression data }
 \item{genotypeFile}{ file containing genotypical data for children }
 \item{correction}{ specifies whether expression data should be corrected (if there are easily distinguishable groups in corelation matrix) }
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
  \code{\link{readExpression}}
  \code{\link{correctChildrenExpression}}
  \code{\link{selectChildrenExpression}}
  \code{\link{readChildrenGenotypes}}
}

\keyword{manip}
