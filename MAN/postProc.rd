\name{postProcessing}
\alias{postProc}


\title{Postprocessing of cross object.}

\description{
  Postprocessing of an object of class cross in sem-automated fashion
}

\usage{
	postProc(cross,n.linkGroups,max.rf.range=c(0.15,0.30),min.lod.range=c(0,3),verbose=FALSE)
	
}

\arguments{
 \item{cross}{ an object of class cross, containing physical or genetic map}
 \item{n.linkGroups}{ expected number of linkage groups}
 \item{max.rf.range}{ range, within which max.rf parameter of formLinkageGroup will be checked}
 \item{min.lod.range}{ range, within which min.lod parameter of formLinkageGroup will be checked}
  \item{verbose}{ be verbose}
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
