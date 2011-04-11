\name{readFiles}
\alias{readFiles}

\title{Reading geno/phenotypic files.}

\description{
  Reads geno/phenotypic files into R environment into special object.
}

\usage{
	readFiles(rils="children",parental="parental",sep="",verbose=FALSE,debugMode=0)
}

\arguments{
 \item{rils}{ Core used to specify names of children phenotypic ("core_phenotypes.txt") and genotypic ("core_genotypes.txt") files.}
 \item{parental}{ Core used to specify names of parental phenotypic ("core_phenotypes.txt") file. }
 \item{sep}{ Separator of values in files. Passed directly to read.table, so "" is a wildcard meaning whitespace.}
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information }
}

\value{
  Object of type ril, containing in this step:
  + essential - $rils$phenotypes - children phenotypic file data
  + optional - $parental$phenotypes - parental phenotypic file data
  + optional - $rils$genotypes$read - children genotypic file data
}

\details{
  To be sure your analysis runs smoothly, data files should have specific format - data separated by tabs, first row provides colnames and first column rownames for the matrix. No addictional remarks, header etc.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
	Under tender patronage of: Danny Arends \email{Danny.Arends@gmail.com}
}

\examples{
	setwd(paste(.Library,"pheno2geno/data",sep="/"))
	ril <- readFiles()
	ril$parental$phenotypes
}

\seealso{
  \code{\link{preprocessData}}
  \code{\link{toGenotypes}}
}

\keyword{manip}
