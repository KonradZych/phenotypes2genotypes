\name{fakePopulation}
\alias{fakePopulation}


\title{Faking population object.}

\description{
  Simulates basic population object for use in examples.
}

\usage{
  fakePopulation(n.founders = 4, n.offspring = 100, n.markers=100,n.chromosomes=10, type = c("riself", "f2", "bc", "risib"), n.mixups=0, verbose=FALSE,...)
}

\arguments{
 \item{n.founders}{ number of founders to be simulated}
 \item{n.offspring}{ number of offspring individuals to be simulated}
 \item{n.markers}{ number of markers individuals to be simulated}
 \item{n.chromosomes}{ number of chromosomes individuals to be simulated}
 \item{type}{ type of the cross to be faked:
   \itemize{
    \item{riself}{ - RILs by selfing}
    \item{f2}{ - f2 cross}
    \item{bc}{ - back cross}
    \item{risib}{ - RILs by siblinmg mating}
  } 
 }
 \item{n.mixups}{ number of mixups to be faked}
 \item{verbose}{ be verbose}
 \item{...}{ to be passed to sim.cross}
}

\value{
  An object of class \code{\link{population}}.
}

\details{
  TODO
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	population <- fakePopulation()
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
  \code{\link{readFiles}}
  \code{\link{intoPopulation}}
  \code{\link{createPopulation}}
  \itemize{
    \item \code{\link[qtl]{sim.cross}} - function from R/qtl package used to simulate genotypic data
    \item \code{\link{readFiles}} - reading files from disc into object of class population
    \item \code{\link{intoPopulation}} - creating object of class population from data already in R environment
    \item \code{\link{createPopulation}} - add data to existing object of class population
  }
}

\keyword{manip}
