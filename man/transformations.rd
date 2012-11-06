\name{transform}
\alias{transform}

\title{Basic functions to do transformation / normalization of phenotypes.}

\description{
  Basic functions to do transformation / normalization of phenotypes.
}

\usage{
  transform(x, transformations=c("nothing","log","sqrt","reciprocal","probit","logit"), ..., verbose=TRUE)
}

\arguments{
  \item{x}{ data matrix with measurements, Rows: Traits/Phenotypes columns: Individuals. }
  \item{transformations}{
    which function should be used to transform the data:
    \itemize{
      \item{nothing}{ - no data transformation performed.}
      \item{\code{\link[base]{log}}}{ - log(data)}
      \item{\code{\link[base]{sqrt}}}{ - sqrt(data)}
      \item{reciprocal}{ - 1/(data)}
      \item{\code{\link[VGAM]{probit}}}{ - probit transformation}
      \item{\code{\link[VGAM]{logit}}}{ -  logit transformation}
    }
  }
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
