\name{Create new map}
\alias{createNewMap}


\title{Creating de novo genetic map.}

\description{
  Creating de novo genetic map.
}

\usage{
	createNewMap(population, orderUsing=c("none","map_genetic","map_physical"), n.chr, use=c("geno","rf"),verbose=FALSE, debugMode=0)
	
}

\arguments{
 \item{population}{ Population type object, must contain parental phenotypic data.}
  \item{orderUsing}{ 
  which map should be used to order markers (Default - none, so markers are all put in 1 chromosome, with distance 1 cM between)
  \itemize{
    \item{map_genetic}{ - supplied genetic map}
    \item{map_physical}{ - supplied physical map}
  }
  }
 \item{n.chr}{ range, within which max.rf parameter of formLinkageGroup will be checked}
 \item{use}{ what kind of data should be used for splitting genotypes into chromosomes
   \itemize{
    \item{geno}{ - genotypes}
    \item{rf}{ - recombination fractions}
  }
 }
 \item{verbose}{ be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
}

\value{
  an object of class cross
}

\details{
postProc function is makign use of most basic piece of inromation possible which is number of linkage groups
(normally equall to number of chromosome) expected. It is using \code{\link{formLinkageGroups}} functions from
R/qtl package with number of different parameter values and afterwards assesses which combination was the best one.
}
\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	#population <- fakePopulation()
	#population <- findDiffExpressed(population)
	#cross <- toGenotypes(population)
	#cross <- postProc(cross,10)
}

\seealso{
  \code{\link{orderChromosomes}} - ordering chromosomes of an object of class cross using majority rule
  \code{\link{rearrangeMarkers}} - rearrangeing markers inside an object of class cross using correlation
  \code{\link{reduceChromosomesNumber}} - removing all but certain number of chromosomes from an object of class cross
}

\keyword{manip}
