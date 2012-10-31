\name{transform}
\alias{transform}

\title{Basic functions to do transformation / normalization of phenotypes.}

\description{
  Basic functions to do transformation / normalization of phenotypes.
}

\usage{
	transform <- function(matrix, transformations=c("nothing","log","sqrt","reciprocal","probit","logit"), ... , verbose=TRUE){
}

\arguments{
  \item{matrix}{ data matrix with measurements, Rows: Traits/Phenotypes columns: Individuals. }
  \item{...}{ Passed to the underlying test function. }
  \item{verbose}{ Be verbose.}
}

\value{
	List with matrices.
}

\author{
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com}
}

\examples{
	data <- matrix(runif(1000),10,100)
  resA <- transform(data, c("log","logit"))
  resB <- transform(data, c("reciprocal","probit"))
}

\seealso{
  \itemize{
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or chromosome assignment vector.}
  }
}

\keyword{manip}
