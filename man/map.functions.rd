\name{map.functions}
\alias{avg_map_distances}
\alias{map_distances}

\title{Functions to provide some descriptive statistics on genetic maps}

\description{
  Functions to provide some descriptive statistics on genetic maps
}

\usage{
  avg_map_distances(m)
  map_distances(m)
}

\arguments{
  \item{m}{ An object of class \code{cross} or \code{map}, See \code{\link[qtl]{read.cross}} or \code{\link[qtl]{pull.map}} for details. }
}

\value{
  A list with per chromosomes either the average map distance or the total distance
}

\author{
	Danny Arends \email{Danny.Arends@gmail.com}, Konrad Zych \email{k.zych@rug.nl}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com}
}

\examples{
	data(yeastCross)
  avg_map_distances(yeastCross)
  map_distances(yeastCross)
}

\seealso{
  \itemize{
    \item{\code{\link{read.population}}}{ - Load genotype, phenotype, genetic map data files into R environment into a population object.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or vector showing how chromosomes should be assigned.}
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{find.diff.expressed}}}{ - Using Rank Product or student t-test analysis to select differentially expressed genes.}
  }
}

\keyword{manip}
