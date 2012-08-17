\name{fake.population}
\alias{fake.population}

\title{Simulate a population object.}

\description{
  Simulates a basic population object for use in examples.
}

\usage{
  fake.population(n.founders = 4, n.offspring = 100, n.markers=100,n.chromosomes=10, type = c("riself", "f2", "bc", "risib"), n.mixups=0, verbose=FALSE,...)
}

\arguments{
 \item{n.founders}{ Number of founders to be simulated.}
 \item{n.offspring}{ Number of offspring individuals to be simulated.}
 \item{n.markers}{ Number of markers individuals to be simulated.}
 \item{n.chromosomes}{ Number of chromosomes individuals to be simulated.}
 \item{type}{ Type of the cross to be faked:
   \itemize{
    \item{riself}{ - RILs by selfing.}
    \item{f2}{ - f2 cross.}
    \item{bc}{ - back cross.}
    \item{risib}{ - RILs by sibling mating.}
  } 
 }
 \item{n.mixups}{ Number of mixups to be faked.}
 \item{verbose}{ Be verbose.}
 \item{...}{ To be passed to \code{\link[qtl]{sim.cross}}.}
}

\value{
  An object of class \code{\link{population}}. See \code{\link{create.population}} for details.
}

\details{
  This function simulates a population object that can be used for further analysis.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	population <- fake.population()
}

\references{
  Copenhaver, G. P., Housworth, E. A. and Stahl, F. W. (2002) Crossover
  interference in arabidopsis.  \emph{Genetics} \bold{160}, 1631--1639.

  Foss, E., Lande, R., Stahl, F. W. and Steinberg, C. M. (1993) Chiasma
  interference as a function of genetic distance. \emph{Genetics}
  \bold{133}, 681--691.

  Zhao, H., Speed, T. P. and McPeek, M. S. (1995) Statistical analysis
  of crossover interference using the chi-square model.  \emph{Genetics}
  \bold{139}, 1045--1056.

  Broman, K. W. (2005) The genomes of recombinant inbred lines
  \emph{Genetics} \bold{169}, 1133--1146.

  Teuscher, F. and Broman, K. W. (2007) Haplotype probabilities for
  multiple-strain recombinant inbred lines.  \emph{Genetics} \bold{175},
  1267--1274. 
}

\seealso{
  \itemize{
    \item{\code{\link{read.population}}}{ - Load genotype, phenotype, genetic map data files into R environment into a population object.}
    \item{\code{\link{add.to.population}}}{ - Add data to existing population object.}
    \item{\code{\link[qtl]{sim.cross}}}{ - Function from R/qtl package used to simulate genotypic data.}
    \item{\code{\link{create.population}}}{ - Create object of class population from data already in R environment.}
  }
}

\keyword{manip}
