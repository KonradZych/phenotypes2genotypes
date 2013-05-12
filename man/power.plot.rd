\name{power.plot }
\alias{power.plot }

\title{Comparison of power of qtl detection.}

\description{
  Plots maximal values of QTL peak measured on the same phenotypes in two crosses.
}

\usage{
  power.plot(cross1,cross2,qtlThr=5,nPheno=500,verbose=FALSE,...)
}

\arguments{
 \item{cross1}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details.}
 \item{cross2}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details.}
 \item{qtlThr}{ Threshold for assessing the significance of the QTL peak.}
 \item{nPheno}{ Nr of phenotypes that will be scanned for QTLs. Phenotypes are selected randomly.}
 \item{verbose}{ Be verbose. }
 \item{...}{ Arguments passed to scanone function (see \code{\link[qtl]{scanone}}).}
}

\details{
Plots maximal values of QTL peak measured on the same phenotypes in two crosses. This give a good comparison of power to detect the QTLs between crosses, if the number of phenotypes scanned is large enough.
}

\value{
	List of maximal values for all the scanned phenotypes.
}

\author{
	Konrad Zych \email{k.zych@rug.nl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{k.zych@rug.nl}
}

\examples{
	#TO ADD
}

\seealso{
  \itemize{
    \item{\code{\link{plotMapComparison}}}{ -  Plotting routine for comparison of two genetic maps.}
    \item{\code{\link{projectOldMarkers}}}{ -   Plotting routine for showing how markers from original map are placed on saturated map.}
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or chromosome assignment vector.}
}
}

\keyword{manip}
