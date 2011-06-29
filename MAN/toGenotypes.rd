\name{toGenotypes}
\alias{toGenotypes}
\alias{use}
\alias{treshold}
\alias{overlapInd}
\alias{proportion}
\alias{margin}
\alias{splitMethod}
\alias{numberOfChromosomes}

\title{Creating genotypes from children phenotypes.}

\description{
  Creating genotypes from children phenotypes using parental data and saving cross object.
}

\usage{
	toGenotypes(population, genotype=c("simulated","real"), orderUsing=c("map_genetic","map_physical"), splitMethod=c("EM","mean"),treshold=0.01, overlapInd = 0, proportion = c(50,50), margin = 15, numberOfChromosomes = NULL, verbose=FALSE, debugMode=0,...)
}

\arguments{
 \item{population}{ Population type object, must contain parental phenotypic data.}
 \item{genotype}{ 
	Which genotypic matrix should be saved to file:
	\itemize{
	\item{simulated}{ - made by toGenotypes}
	\item{real}{ - supported by user and read from file}
	}
}
 \item{orderUsing}{ 
	which map should be used to order markers (by default - none, so markers are all put in 1 chromosome, with distance 1 cM between)
	\itemize{
	\item{map_genetic}{ - simulated data orderd using supported genetic map}
	\item{map_physical}{ - simulated data orderd using supported physical map}
	}
}

 \item{splitMethod}{ Splitting markers using founders mean value or more sofisticated fitting of normal distributions by EM algoritm.}
 \item{treshold}{ If Rank Product pval for gene is lower that this value, we assume it is being diff. expressed.}
 \item{overlapInd}{ Number of individuals that are allowed in the overlap }
 \item{proportion}{ Proportion of individuals expected to carrying a certain genotype }
 \item{margin}{ Proportion is allowed to varry between this margin (2 sided) }
 \item{numberOfChromosomes}{ How many chromosomes should map contain - leaving only that number of linkage groups and dropping one with higher numbers, this can be dangerous, so be sure you know what you're doing. }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print our checks, 2: print additional time information }
  \item{...}{ Parameters passed to \code{\link[qtl]{formLinkageGroups}}. }
}

\value{
  Cross type object. 
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
\dontrun{
	setwd(paste(.Library,"pheno2geno/data",sep="/"))
	ril <- readFiles()
	ril <- preprocessData(ril)
	#cross <- toGenotypes(ril,use="simulated",minChrLength=0,treshold=0.5,margin=50,max.rf=10)
	}
}

\seealso{
  \code{\link{readFiles}}
  \code{\link{preprocessData}}
}

\keyword{manip}
